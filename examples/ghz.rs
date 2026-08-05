//! GHZ state on N qubits: (|0…0> + |1…1>) / sqrt(2).
//!
//! Run with: `cargo run --example ghz`

use qsimulator::Circuit;

fn main() {
    let n = 4;

    let mut c = Circuit::new(n);
    c.h(0);
    for target in 1..n {
        c.cnot(0, target);
    }

    println!("{}\n", c.diagram());

    let counts = c.sample(1000, 7);
    let all_ones = (1usize << n) - 1;

    let mut total = 0;
    for (&outcome, &count) in &counts {
        assert!(
            outcome == 0 || outcome == all_ones,
            "GHZ collapsed to an unexpected outcome {outcome:0n$b}",
            n = n
        );
        total += count;
        println!("|{outcome:0n$b}>: {count}", n = n);
    }
    assert_eq!(total, 1000);
}
