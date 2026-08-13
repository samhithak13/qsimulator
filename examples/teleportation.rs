//! Quantum teleportation. An arbitrary single-qubit state prepared on qubit 0
//! is moved to qubit 2 using a shared Bell pair, a Bell-basis measurement, and
//! classically conditioned corrections.
//!
//! This uses the lower-level `State` API directly. `Circuit` models the
//! mid-circuit measurement (`Circuit::measure`), but not the *feed-forward*:
//! the corrections below depend on which outcomes came up, and a `Circuit`
//! has no classical bits to branch on. Reading the outcomes back needs the
//! `State` API, where `measure_qubit` returns them.
//!
//! Run with: `cargo run --example teleportation`

use qsimulator::{gates, Rng, State};

fn main() {
    // The state to teleport: Ry(theta)|0>, whose |1> probability is sin^2(theta/2).
    let theta = 0.7_f64;
    let expected_p1 = (theta / 2.0).sin().powi(2);

    let mut state = State::new(3);
    state.apply_1q(&gates::ry(theta), 0); // qubit 0 = payload

    // Entangle qubits 1 and 2 into a Bell pair.
    state.apply_1q(&gates::h(), 1);
    state.apply_controlled_1q(&gates::x(), 1, 2);

    // Rotate qubits 0 and 1 into the Bell measurement basis.
    state.apply_controlled_1q(&gates::x(), 0, 1);
    state.apply_1q(&gates::h(), 0);

    // Measure the two control qubits.
    let mut rng = Rng::new(2024);
    let m0 = state.measure_qubit(0, &mut rng);
    let m1 = state.measure_qubit(1, &mut rng);

    // Conditioned corrections recover the payload on qubit 2.
    if m1 {
        state.apply_1q(&gates::x(), 2);
    }
    if m0 {
        state.apply_1q(&gates::z(), 2);
    }

    let p1 = state.prob_qubit_one(2);
    println!("measured (m0, m1) = ({m0}, {m1})");
    println!("qubit 2 P(|1>) = {p1:.6}  (expected {expected_p1:.6})");

    assert!(
        (p1 - expected_p1).abs() < 1e-9,
        "teleported state does not match the original"
    );
    println!("teleportation succeeded");
}
