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

/// General single-qubit gate U3(θ, φ, λ), the OpenQASM `u3`:
/// `[[cos(θ/2),        -e^{iλ}·sin(θ/2)],`
/// ` [e^{iφ}·sin(θ/2),  e^{i(φ+λ)}·cos(θ/2)]]`.
///
/// Every single-qubit unitary is a `u3` up to global phase — e.g.
/// `u3(π,0,π)` = X, and `u3(0,0,λ)` is the phase gate [`p`]`(λ)`.
pub fn u3(theta: f64, phi: f64, lambda: f64) -> Gate {
    let (st, ct) = (theta / 2.0).sin_cos();
    let m00 = c(ct, 0.0);
    // -e^{iλ}·sin(θ/2)
    let m01 = c(-lambda.cos() * st, -lambda.sin() * st);
    // e^{iφ}·sin(θ/2)
    let m10 = c(phi.cos() * st, phi.sin() * st);
    // e^{i(φ+λ)}·cos(θ/2)
    let (spl, cpl) = (phi + lambda).sin_cos();
    let m11 = c(cpl * ct, spl * ct);
    [[m00, m01], [m10, m11]]
}

/// Single-qubit gate U2(φ, λ) = U3(π/2, φ, λ), the OpenQASM `u2`.
///
/// For example `u2(0, π)` = H.
pub fn u2(phi: f64, lambda: f64) -> Gate {
    u3(std::f64::consts::FRAC_PI_2, phi, lambda)
}

/// Decompose a 2x2 unitary `g` into a global phase and U3 angles such that
/// `g == e^{i·γ}·u3(θ, φ, λ)`. Returns `(γ, θ, φ, λ)`.
///
/// Used by OpenQASM export to render an arbitrary controlled-U as a phase on
/// the control (`u1(γ)`) followed by `cu3(θ, φ, λ)`, since the OpenQASM subset
/// has no gate for a controlled arbitrary unitary directly.
pub(crate) fn u3_decompose(g: &Gate) -> (f64, f64, f64, f64) {
    const EPS: f64 = 1e-12;
    // |g00| = cos(θ/2), |g10| = sin(θ/2).
    let cos_half = g[0][0].norm();
    let sin_half = g[1][0].norm();
    let theta = 2.0 * sin_half.atan2(cos_half);

    if cos_half < EPS {
        // θ ≈ π: g00 ≈ 0, so γ is a gauge freedom — fix γ = 0.
        // g10 = e^{i(γ+φ)}, g01 = -e^{i(γ+λ)}.
        (0.0, theta, g[1][0].arg(), (-g[0][1]).arg())
    } else if sin_half < EPS {
        // θ ≈ 0: diagonal, φ is a gauge freedom — fix φ = 0.
        // g00 = e^{iγ}, g11 = e^{i(γ+λ)}.
        let gamma = g[0][0].arg();
        (gamma, theta, 0.0, g[1][1].arg() - gamma)
    } else {
        // Generic: γ from g00, then φ from g10 and λ from -g01.
        let gamma = g[0][0].arg();
        (
            gamma,
            theta,
            g[1][0].arg() - gamma,
            (-g[0][1]).arg() - gamma,
        )
    }
}

/// The square root of X: `sqrt_x() · sqrt_x() == x()`.
///
/// Used by OpenQASM export as the `V` of the Barenco square-root recursion,
/// which is how a multi-controlled X is decomposed when the register has no
/// spare qubit to borrow.
pub(crate) fn sqrt_x() -> Gate {
    [[c(0.5, 0.5), c(0.5, -0.5)], [c(0.5, -0.5), c(0.5, 0.5)]]
}

/// The inverse of [`sqrt_x`]: `sqrt_x_dg() · sqrt_x() == I`.
pub(crate) fn sqrt_x_dg() -> Gate {
    [[c(0.5, -0.5), c(0.5, 0.5)], [c(0.5, 0.5), c(0.5, -0.5)]]
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

    #[test]
    fn u3_specializes_to_known_gates() {
        let id = [[c(1.0, 0.0), c(0.0, 0.0)], [c(0.0, 0.0), c(1.0, 0.0)]];
        assert_gate_eq(&u3(0.0, 0.0, 0.0), &id);
        assert_gate_eq(&u3(PI, 0.0, PI), &x());
        assert_gate_eq(&u3(0.0, 0.0, 0.9), &p(0.9));
    }

    #[test]
    fn u2_zero_pi_is_hadamard() {
        assert_gate_eq(&u2(0.0, PI), &h());
    }

    /// Reconstruct `e^{iγ}·u3(θ, φ, λ)`.
    fn recompose(gamma: f64, theta: f64, phi: f64, lambda: f64) -> Gate {
        let phase = c(gamma.cos(), gamma.sin());
        let g = u3(theta, phi, lambda);
        [
            [phase * g[0][0], phase * g[0][1]],
            [phase * g[1][0], phase * g[1][1]],
        ]
    }

    #[test]
    fn u3_decompose_reconstructs_every_gate() {
        // Named gates, rotations, diagonal (θ=0) and anti-diagonal (θ=π) cases,
        // and a few generic unitaries must all round-trip through the
        // decomposition.
        let cases = [
            x(),
            y(),
            z(),
            h(),
            s(),
            t(),
            sdg(),
            tdg(),
            rz(0.7),
            p(1.3),
            u3(0.5, -0.6, 0.7),
            u3(PI, 0.4, -0.9),                                        // θ = π
            u3(0.0, 0.2, 1.1),                                        // θ = 0 (diagonal)
            [[c(1.0, 0.0), c(0.0, 0.0)], [c(0.0, 0.0), c(1.0, 0.0)]], // identity
        ];
        for g in cases {
            let (gamma, theta, phi, lambda) = u3_decompose(&g);
            assert_gate_eq(&recompose(gamma, theta, phi, lambda), &g);
        }
    }

    #[test]
    fn sqrt_x_squares_to_x_and_inverts() {
        let id = [[c(1.0, 0.0), c(0.0, 0.0)], [c(0.0, 0.0), c(1.0, 0.0)]];
        assert_gate_eq(&mul(&sqrt_x(), &sqrt_x()), &x());
        assert_gate_eq(&mul(&sqrt_x_dg(), &sqrt_x()), &id);
    }

    #[test]
    fn u3_is_unitary() {
        let g = u3(0.7, -1.1, 2.3);
        let dot = |i: usize, j: usize| g[i][0] * g[j][0].conj() + g[i][1] * g[j][1].conj();
        let prod = [[dot(0, 0), dot(0, 1)], [dot(1, 0), dot(1, 1)]];
        let id = [[c(1.0, 0.0), c(0.0, 0.0)], [c(0.0, 0.0), c(1.0, 0.0)]];
        assert_gate_eq(&prod, &id);
    }
}
