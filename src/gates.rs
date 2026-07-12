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
