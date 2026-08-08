//! Quantum Fourier transform on N qubits.
//!
//! The QFT maps the all-zeros state to the uniform superposition — every basis
//! state with equal, real, positive amplitude `1/sqrt(2^n)` — which this
//! example builds and checks.
//!
//! Run with: `cargo run --example qft`

use qsimulator::Circuit;

/// Append an n-qubit QFT to `c`: a Hadamard on each qubit, controlled phase
/// rotations between every pair, and a final bit reversal.
fn qft(c: &mut Circuit, n: usize) {
    for j in 0..n {
        c.h(j);
        for k in (j + 1)..n {
            // Controlled phase of 2*pi / 2^(k-j+1) from qubit k onto qubit j.
            let angle = std::f64::consts::PI / (1u64 << (k - j)) as f64;
            c.cp(angle, k, j);
        }
    }
    // Reverse the qubit order.
    for i in 0..n / 2 {
        c.swap(i, n - 1 - i);
    }
}

fn main() {
    let n = 4;
    let mut c = Circuit::new(n);
    qft(&mut c, n); // applied to |0...0>

    let state = c.run();
    println!("{}\n", c.diagram());

    let dim = 1usize << n;
    let expected = 1.0 / (dim as f64).sqrt();
    for i in 0..dim {
        let amp = state.amplitudes()[i];
        // QFT|0> is the uniform superposition: every amplitude is real 1/sqrt(N).
        assert!(
            (amp.re - expected).abs() < 1e-9 && amp.im.abs() < 1e-9,
            "amplitude {i} = {amp:?}, expected {expected}"
        );
    }
    println!("QFT|0...0> is the uniform superposition: all {dim} amplitudes = {expected:.4}");
}
