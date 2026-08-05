//! Grover's search on two qubits. With a four-element search space, a single
//! Grover iteration rotates the uniform superposition exactly onto the marked
//! state, so it is found with certainty.
//!
//! Run with: `cargo run --example grover`

use qsimulator::Circuit;

fn main() {
    let mut c = Circuit::new(2);

    // Uniform superposition over the four basis states.
    c.h(0).h(1);

    // Oracle: flip the phase of the marked state |11>.
    c.cz(0, 1);

    // Diffusion operator (inversion about the mean).
    c.h(0).h(1).x(0).x(1).cz(0, 1).x(0).x(1).h(0).h(1);

    let state = c.run();
    println!("{}\n", c.diagram());
    for i in 0..4 {
        println!("|{i:02b}>: p = {:.3}", state.probability(i));
    }

    let marked = state.probability(0b11);
    assert!(
        (marked - 1.0).abs() < 1e-9,
        "expected the marked state with probability 1, got {marked}"
    );
    println!("\nfound |11> with probability {marked:.3}");
}
