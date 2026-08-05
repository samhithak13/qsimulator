# Benchmarks

`throughput.rs` measures single-gate application throughput — the inner loop
of the simulator. It is a dependency-free custom harness (`harness = false`).

## Running

```bash
cargo bench
```

Each case applies one gate to a full `n`-qubit register, auto-scales its
iteration count until it has run for at least ~300 ms, and reports the mean
time per call and the derived amplitude update rate (`2^n / time`).

## Protocol

- Release build (`cargo bench` compiles with optimizations).
- One warm-up call, then geometric iteration growth to the time floor.
- Wall-clock timing via `std::time::Instant`; results are passed through
  `std::hint::black_box` so the work is not optimized away.

Numbers are machine- and load-dependent. Treat them as relative, not
absolute, and regenerate them on the target machine.

## Representative results

Measured on an Apple M1 (aarch64), release build. Rates in millions of
amplitude updates per second (higher is better).

| Case (n = 22)             | Before | After |
|---------------------------|-------:|------:|
| `apply_1q`, target = 0    |  ~201  | ~320  |
| `apply_1q`, target = n-1  |  ~240  | ~325  |
| `apply_controlled_1q`     |  ~338  | ~415  |

"Before" is the original index-based kernel; "After" restructures `apply_1q`
into a bounds-check-free walk over `chunks_exact_mut` split into the target
bit's `|0>`/`|1>` halves, and hoists the gate-matrix entries into locals in
all three kernels. The largest gain is on a low target qubit, where the old
kernel's small stride blocks had the most per-iteration overhead.
