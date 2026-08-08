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

## Parallel kernels

Building with `--features parallel` runs the single-qubit and (multi-)
controlled kernels across threads (rayon). Same machine, `n = 22`:

```bash
cargo bench --features parallel
```

| Case (n = 22)            | Serial | Parallel |
|--------------------------|-------:|---------:|
| `apply_1q`, target = 0   |  ~500  |   ~760   |
| `apply_1q`, target = n-1 |  ~510  |  ~1200   |
| `apply_controlled_1q`    |  ~820  |  ~1000   |

Each kernel parallelizes over blocks for a low target qubit and within a block
for a high one. The speedup is sublinear in the eight cores because a
state-vector sweep is memory-bandwidth bound; the high-target single-qubit
case, which splits into two large contiguous halves, scales best (~2.4x). The
controlled kernel gains less (~1.2x): it is branchy and touches only the pairs
whose control bits are set.
