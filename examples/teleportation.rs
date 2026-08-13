//! Quantum teleportation. An arbitrary single-qubit state prepared on qubit 0
//! is moved to qubit 2 using a shared Bell pair, a Bell-basis measurement, and
//! corrections conditioned on what that measurement found.
//!
//! Written with the `Circuit` builder, which models the whole protocol:
//! `measure` collapses and records an outcome, and `if_classical_eq` guards the
//! corrections on it. The interesting property is that teleportation works for
//! *every* measurement outcome — the corrections are exactly what undoes each
//! one — so this runs several seeds and checks all of them.
//!
//! Run with: `cargo run --example teleportation`

use qsimulator::Circuit;

/// Build the protocol for a payload state `Ry(theta)|0>` on qubit 0.
fn teleport(theta: f64) -> Circuit {
    let mut c = Circuit::new(3);

    // Qubit 0 holds the payload; qubits 1 and 2 share a Bell pair.
    c.ry(theta, 0);
    c.h(1).cnot(1, 2);

    // Rotate qubits 0 and 1 into the Bell measurement basis and read them out.
    // Qubit i is recorded in classical bit i, so the register holds
    // m0 + 2*m1 once both have been measured.
    c.cnot(0, 1).h(0);
    c.measure(0).measure(1);

    // Undo whichever Bell outcome came up. m1 calls for an X on the target and
    // m0 for a Z, so the four register values need: nothing, Z, X, then both.
    c.if_classical_eq(0b10, |b| {
        b.x(2);
    });
    c.if_classical_eq(0b01, |b| {
        b.z(2);
    });
    c.if_classical_eq(0b11, |b| {
        b.x(2).z(2);
    });
    c
}

fn main() {
    let theta = 0.7_f64;
    let expected_p1 = (theta / 2.0).sin().powi(2);
    let circuit = teleport(theta);

    println!("{circuit}\n");
    println!("payload P(|1>) = {expected_p1:.6}");

    // Different seeds take different branches of the Bell measurement. Walk
    // seeds until all four have come up, checking each one: a correction that
    // was wrong for a single branch would otherwise hide behind the three that
    // happened to be sampled.
    let mut seen = [false; 4];
    for seed in 0..64 {
        let state = circuit.run_seeded(seed);
        let p1 = state.prob_qubit_one(2);
        // Which branch this seed landed on, read off the collapsed controls.
        let m0 = state.prob_qubit_one(0) > 0.5;
        let m1 = state.prob_qubit_one(1) > 0.5;
        let branch = usize::from(m0) | usize::from(m1) << 1;

        assert!(
            (p1 - expected_p1).abs() < 1e-9,
            "teleported state wrong on branch ({m0}, {m1}): {p1} vs {expected_p1}"
        );
        if !seen[branch] {
            seen[branch] = true;
            println!("branch ({m0}, {m1}) via seed {seed}: qubit 2 P(|1>) = {p1:.6}");
        }
        if seen.iter().all(|s| *s) {
            break;
        }
    }
    assert!(seen.iter().all(|s| *s), "not every branch was exercised");

    println!("\nteleportation succeeded on all four branches");
}
