//! Noise channels: how a circuit degrades, and how fast.
//!
//! Channels are simulated by quantum trajectories — each shot samples one
//! Kraus operator — so a single run is one sample and the numbers below come
//! from `sample`. That also means they carry sampling error, which is why the
//! tolerances are loose rather than exact.
//!
//! Run with: `cargo run --example noise`

use qsimulator::Circuit;

const SHOTS: usize = 20_000;
const SEED: u64 = 42;

/// Fraction of shots in which qubit `q` read |1>.
fn fraction_one(c: &Circuit, q: usize) -> f64 {
    let hist = c.sample(SHOTS, SEED);
    let ones: usize = hist
        .iter()
        .filter(|(state, _)| *state >> q & 1 == 1)
        .map(|(_, n)| n)
        .sum();
    ones as f64 / SHOTS as f64
}

/// Fraction of shots where the two qubits of a Bell pair disagreed.
fn disagreement(c: &Circuit) -> f64 {
    let hist = c.sample(SHOTS, SEED);
    let split: usize = hist
        .iter()
        .filter(|(state, _)| (*state & 1) != (*state >> 1 & 1))
        .map(|(_, n)| n)
        .sum();
    split as f64 / SHOTS as f64
}

fn main() {
    // A noiseless Bell pair always agrees. Depolarizing noise on one half
    // breaks that: the channel applies X, Y or Z each with p/3, and the two
    // that flip the bit (X and Y) show up as disagreement, so the rate climbs
    // towards 2p/3.
    println!("Bell pair, depolarizing noise on qubit 1");
    println!("    p     disagree   expected (2p/3)");
    for p in [0.0, 0.05, 0.2, 0.5, 0.9] {
        let mut c = Circuit::new(2);
        c.h(0).cnot(0, 1).depolarizing(p, 1);
        let measured = disagreement(&c);
        let expected = 2.0 * p / 3.0;
        println!("  {p:4.2}     {measured:.4}     {expected:.4}");
        assert!(
            (measured - expected).abs() < 0.02,
            "depolarizing({p}): {measured} vs {expected}"
        );
    }

    // Amplitude damping is T1 relaxation: an excited qubit decays to |0> with
    // probability gamma each time it is applied, so surviving k rounds has
    // probability (1-gamma)^k — an exponential, the shape a T1 measurement
    // actually produces on hardware.
    let gamma = 0.2;
    println!("\nAmplitude damping, gamma = {gamma}: survival of |1> over rounds");
    println!("  rounds   survived   expected ((1-g)^k)");
    for rounds in 0..6 {
        let mut c = Circuit::new(1);
        c.x(0);
        for _ in 0..rounds {
            c.amplitude_damping(gamma, 0);
        }
        let measured = fraction_one(&c, 0);
        let expected = (1.0 - gamma).powi(rounds);
        println!("  {rounds:6}     {measured:.4}     {expected:.4}");
        assert!(
            (measured - expected).abs() < 0.02,
            "{rounds} rounds: {measured} vs {expected}"
        );
    }

    // Phase damping takes no energy, so populations never move — but it
    // destroys the interference that makes H twice an identity.
    println!("\nPhase damping between two Hadamards (no damping returns |0>)");
    println!("  gamma    P(|1>)");
    for gamma in [0.0, 0.5, 1.0] {
        let mut c = Circuit::new(1);
        c.h(0).phase_damping(gamma, 0).h(0);
        let measured = fraction_one(&c, 0);
        // Coherence survives as sqrt(1-gamma), so P(|1>) = (1 - sqrt(1-g))/2.
        let expected = (1.0 - (1.0 - gamma).sqrt()) / 2.0;
        println!("  {gamma:5.2}    {measured:.4}   (expected {expected:.4})");
        assert!(
            (measured - expected).abs() < 0.02,
            "phase_damping({gamma}): {measured} vs {expected}"
        );
    }

    println!("\nevery channel matched its analytic rate");
}
