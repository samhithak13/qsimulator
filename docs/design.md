# Design notes

## State-vector model

An `n`-qubit register is represented by a vector of `2^n` complex
amplitudes. Basis states are indexed in **little-endian** order: bit `q`
of the index corresponds to qubit `q`. This makes single-qubit gate
application a stride-`2^q` butterfly over amplitude pairs.

## Gate application

- **Single-qubit gates** apply a 2x2 unitary to every pair of amplitudes
  that differ only in the target bit (`state.apply_1q`).
- **Controlled gates** apply the 2x2 unitary only to pairs where the
  control bit is set (`state.apply_controlled_1q`).
- **SWAP** (`state.swap_qubits`) exchanges two qubits by swapping the
  amplitudes of basis states that differ only in those two bits.
- **Multi-controlled gates** (`state.apply_multi_controlled_1q`) apply a
  2x2 unitary only where *every* control bit is set. Zero controls give an
  unconditional gate, one control matches `apply_controlled_1q`, and two
  controls with X give a Toffoli (`circuit.toffoli`).

Both are done in place, so memory is a single `2^n` vector.

The gate set covers the Paulis (X, Y, Z), Hadamard, the phase gates S and
T, and the continuous rotations `rx(θ)`, `ry(θ)`, `rz(θ)`. Each rotation is
the standard `exp(-i·θ/2·P)` for its Pauli `P`, so e.g. `rx(π)` equals X up
to the global phase `-i`.

## Measurement

Sampling a computational-basis outcome uses the Born rule:
`p(i) = |amplitude(i)|^2`. Post-measurement collapse renormalizes the
surviving amplitudes.

- `State::prob_qubit_one(q)` — probability qubit `q` reads |1> (read-only).
- `State::measure_qubit(q, rng)` — sample a single qubit, collapse the
  register onto the measured subspace, and renormalize.
- `State::measure_all(rng)` — Born-rule sample a full basis state and
  collapse onto it, returning the little-endian index.
- `Circuit::sample(shots, seed)` — run once, then measure independent
  clones of the final state `shots` times into a histogram.

Randomness comes from a seedable, dependency-free `xorshift64` generator
(`rng::Rng`). The seed is mixed through the SplitMix64 finalizer so that
even a zero seed produces a healthy, non-degenerate stream. Because the
generator is fully deterministic in the seed, sampling runs and tests are
exactly reproducible.

## Testing strategy

Every gate has a known truth table / matrix; integration tests assert
exact probabilities for canonical circuits (Bell state, GHZ, single-gate
flips) within a tight epsilon.

---

## Current status & handoff (continue from here)

This section is a running log so any session — including one picked up from
the Claude mobile/web app against this GitHub repo — can continue without
extra context. Update it as features land.

### What works today

Everything below is implemented, tested, and pushed to `main`. CI (fmt +
clippy `-D warnings` + build + test) is green. **43 tests** across 9 test
binaries.

| Area | API | Status |
|---|---|---|
| State vector | `State::new/n_qubits/amplitudes/probability/norm` | ✅ |
| Single-qubit apply | `State::apply_1q` | ✅ |
| Controlled apply | `State::apply_controlled_1q` | ✅ |
| Multi-controlled apply | `State::apply_multi_controlled_1q(gate, controls, target)` | ✅ |
| SWAP | `State::swap_qubits(a, b)` | ✅ |
| Gates | `gates::{x,y,z,h,s,t,rx,ry,rz}` | ✅ |
| Measurement | `State::{prob_qubit_one, measure_qubit, measure_all}` | ✅ |
| Sampling | `Circuit::sample(shots, seed) -> HashMap<usize, usize>` | ✅ |
| RNG | `rng::Rng` (seedable xorshift64, SplitMix64 seeding) | ✅ |
| Circuit builders | `h, x, y, z, s, t, rx, ry, rz, cnot, cz, cu, swap, toffoli, mcx, mcu` | ✅ |
| Circuit display | `Display for Circuit` (ASCII wire diagram) | ✅ |

### File map

- `src/state.rs` — state vector + all `apply_*`, measurement, `swap_qubits`.
- `src/gates.rs` — 2x2 gate matrices (`type Gate = [[Complex64; 2]; 2]`),
  plus a unit-test module for the rotations.
- `src/circuit.rs` — `Op` enum (`Single`, `Controlled`, `Swap`,
  `MultiControlled`), the `Circuit` builder, `run`, and `sample`.
- `src/rng.rs` — the RNG and its unit tests.
- `src/lib.rs` — module wiring and re-exports (`Circuit`, `State`, `Rng`).
- `src/main.rs` — demo: builds a Bell state and samples 1000 shots.
- `tests/` — `bell_state.rs`, `builders.rs`, `ghz.rs`, `measurement.rs`,
  `rotations.rs`, `swap.rs`, `toffoli.rs`.

### Frozen conventions (do not change silently)

- **Little-endian qubit order**: bit `q` of a basis index is qubit `q`.
  A 2-qubit index `0b10` means qubit 1 is |1>, qubit 0 is |0>.
- **Gate storage** is row-major `[[m00, m01], [m10, m11]]`, `Complex64`.
- **Rotations** are `exp(-i·θ/2·P)`; `R_axis(π)` equals its Pauli up to the
  global phase `-i` (tests assert exactly this).
- **RNG is deterministic in its seed** — reproducibility depends on it, so
  keep the xorshift64 recurrence and SplitMix64 seeding stable. Changing
  them will break `same_seed_is_reproducible` and any golden values.
- Controls must differ from the target (asserted); pairs are always
  processed once via the target-bit-0 side.

### Dev commands

```bash
cargo fmt --all -- --check
cargo clippy --all-targets -- -D warnings
cargo build
cargo test
cargo run            # the Bell-state sampling demo
```

New functionality must land with its test in the same commit, and all four
checks above must pass before committing.

### Completed next steps (v0.2)

1. ~~**More builder coverage**~~ — ✅ `y`, `s`, `t`, `cz`, `cu` added with
   tests (12 new tests in `tests/builders.rs`).
2. ~~**Arbitrary controlled-U builder**~~ — ✅ `mcx(controls, target)` and
   `mcu(gate, controls, target)` added with tests (0, 1, and 3 controls).
3. ~~**Circuit inspection / diagram printing**~~ — ✅ `Display for Circuit`
   prints an ASCII wire diagram with gate labels, control dots (`●`), swap
   markers (`×`), and vertical links (`│`). Demo output in `main.rs`.
5. ~~**GHZ + multi-qubit fixtures**~~ — ✅ 3-qubit and 4-qubit GHZ tests
   in `tests/ghz.rs` (probabilities + sampling).

### Suggested next steps (v0.3, not yet started)

Roughly in priority order; none of these exist yet:

4. **Richer CLI** (`main.rs`) — accept a gate list / simple program instead
   of the hard-coded Bell demo.
6. **Performance** — in-place kernels already; consider benchmarks
   and a sparse fast path only after the above.

When adding a gate: add the matrix in `gates.rs`, a builder in `circuit.rs`
(+ `Op` variant and `run` arm if it needs new state machinery), and both a
unit test (matrix/truth-table) and an integration test (circuit behavior).
