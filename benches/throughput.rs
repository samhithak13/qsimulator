//! Gate-application throughput benchmarks.
//!
//! Dependency-free, custom harness (`harness = false`). Each case auto-scales
//! its iteration count until it has run for at least ~300 ms, then reports the
//! mean time per full-register gate application and the derived amplitude
//! update rate.
//!
//! Run with: `cargo bench`
//!
//! Protocol: release build, one warm-up call, geometric iteration growth to
//! the time floor, wall-clock timing via `Instant`. Numbers are machine- and
//! load-dependent; treat them as relative, not absolute.

use std::hint::black_box;
use std::time::{Duration, Instant};

use qsimulator::{gates, State};

const MIN_TIME: Duration = Duration::from_millis(300);

/// Time `f` per call (auto-scaled) and print a row: mean time and, given the
/// register size `n`, the amplitude update rate.
fn bench(name: &str, n: usize, mut f: impl FnMut()) {
    f(); // warm up

    let mut iters: u64 = 1;
    let per = loop {
        let start = Instant::now();
        for _ in 0..iters {
            f();
        }
        let elapsed = start.elapsed();
        if elapsed >= MIN_TIME || iters >= 1 << 30 {
            break elapsed.as_nanos() as f64 / iters as f64;
        }
        iters *= 2;
    };

    let amps = (1u64 << n) as f64;
    let updates_per_sec = amps / (per * 1e-9);
    println!(
        "{name:<26} n={n:<2}  {per:>10.1} ns/call  {:>7.1} M amp/s",
        updates_per_sec / 1e6
    );
}

fn main() {
    let h = gates::h();

    for &n in &[18usize, 20, 22] {
        // Single-qubit gate on the lowest qubit (many small stride blocks).
        let mut state = State::new(n);
        bench("apply_1q H, target=0", n, || state.apply_1q(&h, 0));
        black_box(&state);

        // Single-qubit gate on the highest qubit (one large stride block).
        let mut state = State::new(n);
        let top = n - 1;
        bench("apply_1q H, target=n-1", n, || state.apply_1q(&h, top));
        black_box(&state);

        // A controlled gate across the full register.
        let mut state = State::new(n);
        bench("apply_controlled_1q X", n, || {
            state.apply_controlled_1q(&h, 0, 1)
        });
        black_box(&state);

        println!();
    }
}
