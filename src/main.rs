//! Demo: prepare a Bell state (|00> + |11|)/sqrt(2) and print probabilities.

use qsimulator::Circuit;

fn main() {
    let mut circuit = Circuit::new(2);
    circuit.h(0).cnot(0, 1);
    let state = circuit.run();

    println!("Bell state amplitudes:");
    for (i, amp) in state.amplitudes().iter().enumerate() {
        println!(
            "  |{:02b}>: amp = {:+.3}{:+.3}i, p = {:.3}",
            i,
            amp.re,
            amp.im,
            state.probability(i)
        );
    }
    println!("total probability = {:.6}", state.norm());
}
