//! Standard single-qubit gate matrices.

use num_complex::Complex64;

type Gate = [[Complex64; 2]; 2];

#[inline]
fn c(re: f64, im: f64) -> Complex64 {
    Complex64::new(re, im)
}

/// Pauli-X (NOT) gate.
pub fn x() -> Gate {
    [[c(0.0, 0.0), c(1.0, 0.0)], [c(1.0, 0.0), c(0.0, 0.0)]]
}

/// Pauli-Y gate.
pub fn y() -> Gate {
    [[c(0.0, 0.0), c(0.0, -1.0)], [c(0.0, 1.0), c(0.0, 0.0)]]
}

/// Pauli-Z gate.
pub fn z() -> Gate {
    [[c(1.0, 0.0), c(0.0, 0.0)], [c(0.0, 0.0), c(-1.0, 0.0)]]
}

/// Hadamard gate.
pub fn h() -> Gate {
    let s = std::f64::consts::FRAC_1_SQRT_2;
    [[c(s, 0.0), c(s, 0.0)], [c(s, 0.0), c(-s, 0.0)]]
}

/// Phase gate S = diag(1, i).
pub fn s() -> Gate {
    [[c(1.0, 0.0), c(0.0, 0.0)], [c(0.0, 0.0), c(0.0, 1.0)]]
}

/// T gate = diag(1, e^{i pi/4}).
pub fn t() -> Gate {
    let phase = std::f64::consts::FRAC_PI_4;
    [
        [c(1.0, 0.0), c(0.0, 0.0)],
        [c(0.0, 0.0), c(phase.cos(), phase.sin())],
    ]
}

/// S-dagger gate = S† = diag(1, -i).
pub fn sdg() -> Gate {
    [[c(1.0, 0.0), c(0.0, 0.0)], [c(0.0, 0.0), c(0.0, -1.0)]]
}

/// T-dagger gate = T† = diag(1, e^{-i pi/4}).
pub fn tdg() -> Gate {
    let phase = -std::f64::consts::FRAC_PI_4;
    [
        [c(1.0, 0.0), c(0.0, 0.0)],
        [c(0.0, 0.0), c(phase.cos(), phase.sin())],
    ]
}

/// Phase gate P(λ) = diag(1, e^{iλ}).
///
/// The continuous generalization of the diagonal phase gates: `p(π)` = Z,
/// `p(π/2)` = S, `p(π/4)` = T. Equivalent to OpenQASM's `u1(λ)`.
pub fn p(lambda: f64) -> Gate {
    let (s, co) = lambda.sin_cos();
    [[c(1.0, 0.0), c(0.0, 0.0)], [c(0.0, 0.0), c(co, s)]]
}

/// Rotation about the X axis by angle `theta`:
/// `[[cos(θ/2), -i·sin(θ/2)], [-i·sin(θ/2), cos(θ/2)]]`.
///
/// `rx(π)` equals X up to the global phase `-i`.
pub fn rx(theta: f64) -> Gate {
    let (s, co) = (theta / 2.0).sin_cos();
    [[c(co, 0.0), c(0.0, -s)], [c(0.0, -s), c(co, 0.0)]]
}

/// Rotation about the Y axis by angle `theta`:
/// `[[cos(θ/2), -sin(θ/2)], [sin(θ/2), cos(θ/2)]]`.
///
/// `ry(π)` equals Y up to the global phase `-i`.
pub fn ry(theta: f64) -> Gate {
    let (s, co) = (theta / 2.0).sin_cos();
    [[c(co, 0.0), c(-s, 0.0)], [c(s, 0.0), c(co, 0.0)]]
}

/// Rotation about the Z axis by angle `theta`:
/// `diag(e^{-iθ/2}, e^{iθ/2})`.
///
/// `rz(π)` equals Z up to the global phase `-i`.
pub fn rz(theta: f64) -> Gate {
    let (s, co) = (theta / 2.0).sin_cos();
    [[c(co, -s), c(0.0, 0.0)], [c(0.0, 0.0), c(co, s)]]
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    /// Assert two 2x2 gate matrices are elementwise equal within `1e-12`.
    fn assert_gate_eq(a: &Gate, b: &Gate) {
        for i in 0..2 {
            for j in 0..2 {
                assert!(
                    (a[i][j] - b[i][j]).norm() < 1e-12,
                    "entry ({i},{j}) differs: {:?} vs {:?}",
                    a[i][j],
                    b[i][j]
                );
            }
        }
    }

    /// Multiply `-i` through every entry of `g` (the global phase relating
    /// each `R_axis(π)` to its corresponding Pauli).
    fn neg_i_times(g: &Gate) -> Gate {
        let ni = c(0.0, -1.0);
        [[ni * g[0][0], ni * g[0][1]], [ni * g[1][0], ni * g[1][1]]]
    }

    #[test]
    fn rotations_at_zero_are_identity() {
        let id = [[c(1.0, 0.0), c(0.0, 0.0)], [c(0.0, 0.0), c(1.0, 0.0)]];
        assert_gate_eq(&rx(0.0), &id);
        assert_gate_eq(&ry(0.0), &id);
        assert_gate_eq(&rz(0.0), &id);
    }

    #[test]
    fn rx_pi_is_x_up_to_global_phase() {
        // Rx(π) = -i·X.
        assert_gate_eq(&rx(PI), &neg_i_times(&x()));
    }

    #[test]
    fn ry_pi_is_y_up_to_global_phase() {
        assert_gate_eq(&ry(PI), &neg_i_times(&y()));
    }

    #[test]
    fn rz_pi_is_z_up_to_global_phase() {
        assert_gate_eq(&rz(PI), &neg_i_times(&z()));
    }

    #[test]
    fn rotations_are_unitary() {
        for g in [rx(0.7), ry(-1.3), rz(2.1)] {
            // g · g† should be the identity. Row i of g dotted with the
            // conjugate of row j gives entry (i, j) of the product.
            let dot = |i: usize, j: usize| g[i][0] * g[j][0].conj() + g[i][1] * g[j][1].conj();
            let prod = [[dot(0, 0), dot(0, 1)], [dot(1, 0), dot(1, 1)]];
            let id = [[c(1.0, 0.0), c(0.0, 0.0)], [c(0.0, 0.0), c(1.0, 0.0)]];
            assert_gate_eq(&prod, &id);
        }
    }

    /// Ordinary 2x2 matrix product `a·b`.
    fn mul(a: &Gate, b: &Gate) -> Gate {
        let e = |i: usize, j: usize| a[i][0] * b[0][j] + a[i][1] * b[1][j];
        [[e(0, 0), e(0, 1)], [e(1, 0), e(1, 1)]]
    }

    #[test]
    fn daggers_invert_their_gates() {
        let id = [[c(1.0, 0.0), c(0.0, 0.0)], [c(0.0, 0.0), c(1.0, 0.0)]];
        assert_gate_eq(&mul(&s(), &sdg()), &id);
        assert_gate_eq(&mul(&t(), &tdg()), &id);
    }

    #[test]
    fn phase_gate_specializes_to_z_s_t() {
        assert_gate_eq(&p(PI), &z());
        assert_gate_eq(&p(PI / 2.0), &s());
        assert_gate_eq(&p(PI / 4.0), &t());
    }
}
