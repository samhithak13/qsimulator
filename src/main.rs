//! Demo: prepare a Bell state (|00> + |11>)/sqrt(2), print its amplitudes,
//! then sample it many times to show the measurement statistics.

use qsimulator::Circuit;

fn main() {
    let mut circuit = Circuit::new(2);
    circuit.h(0).cnot(0, 1);

    println!("Circuit:");
    print!("{circuit}");
    println!();

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

    let shots = 1000;
    let seed = 0xC0FF_EE00;
    let histogram = circuit.sample(shots, seed);

    println!("\nSampling {shots} shots (seed = {seed:#x}):");
    for outcome in 0..(1usize << state.n_qubits()) {
        let count = histogram.get(&outcome).copied().unwrap_or(0);
        println!("  |{outcome:02b}>: {count} shots");
    }
}
