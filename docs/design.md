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
clippy `-D warnings` + build + test) is green. **102 tests** (100 unit +
integration across 13 binaries, plus 2 doctests).

| Area | API | Status |
|---|---|---|
| State vector | `State::new/n_qubits/amplitudes/probability/norm` | ✅ |
| Single-qubit apply | `State::apply_1q` | ✅ |
| Controlled apply | `State::apply_controlled_1q` | ✅ |
| Multi-controlled apply | `State::apply_multi_controlled_1q(gate, controls, target)` | ✅ |
| SWAP | `State::swap_qubits(a, b)` | ✅ |
| Gates | `gates::{x,y,z,h,s,t,sdg,tdg,p,u2,u3,rx,ry,rz}` | ✅ |
| Measurement | `State::{prob_qubit_one, measure_qubit, measure_all}` | ✅ |
| Sampling | `Circuit::sample(shots, seed) -> HashMap<usize, usize>` | ✅ |
| RNG | `rng::Rng` (seedable xorshift64, SplitMix64 seeding) | ✅ |
| Circuit single-qubit builders | `h,x,y,z,s,t,sdg,tdg,p,u2,u3,rx,ry,rz` | ✅ |
| Circuit controlled builders | `cnot, cz, crz, cp, cu(gate,c,t), mcx(controls,t), mcu(gate,controls,t)` | ✅ |
| Circuit other builders | `swap, toffoli` (toffoli now delegates to `mcx`) | ✅ |
| Circuit rendering | `Circuit::diagram() -> String` + `Display` (ASCII diagram) | ✅ |
| Text program format | `program::parse(&str) -> Result<Program, String>` | ✅ |
| OpenQASM 2.0 import | `qasm::parse(&str) -> Result<Circuit, String>` (subset) | ✅ |
| OpenQASM 2.0 export | `Circuit::to_qasm() -> Result<String, String>` | ✅ |
| CLI | `qsimulator [FILE\|-\|--emit-qasm FILE\|--help]` | ✅ |

### File map

- `src/state.rs` — state vector + all `apply_*`, measurement, `swap_qubits`.
- `src/gates.rs` — 2x2 gate matrices (`type Gate = [[Complex64; 2]; 2]`),
  plus a unit-test module for the rotations.
- `src/circuit.rs` — `Op` enum (`Single`, `Controlled`, `Swap`,
  `MultiControlled`), the `Circuit` builder, `run`, and `sample`.
- `src/rng.rs` — the RNG and its unit tests.
- `src/program.rs` — text program parser (`parse`, `Program`, `SampleSpec`);
  `parse_angle` is `pub(crate)` and shared with the QASM importer.
- `src/qasm.rs` — OpenQASM 2.0 subset importer (`parse -> Circuit`).
- `src/lib.rs` — module wiring and re-exports (`Circuit`, `State`, `Rng`).
- `src/main.rs` — the CLI: built-in demo, or parse+run a program/`.qasm`/stdin
  (dispatches to `qasm` by `.qasm` extension or `OPENQASM` header).
- `programs/ghz.qsim` — sample text program; `programs/bell.qasm` — sample QASM.
- `tests/` — `bell_state.rs`, `measurement.rs`, `rotations.rs`, `swap.rs`,
  `toffoli.rs`, `single_qubit_builders.rs` (y/s/t), `controlled_builders.rs`
  (cz/cu/mcx/mcu), `diagram.rs` (ASCII rendering), `program.rs` (parser),
  `ghz.rs` (3- and 4-qubit GHZ probabilities + sampling), `qasm.rs`
  (OpenQASM import, oracle-checked against the builder API), `qasm_export.rs`
  (export + round-trip through import).

Note: `Op` variants now carry a `label: &'static str` used only by
`diagram()`; it never affects execution (`run` ignores it). New gate
builders should pass a short static label (e.g. `"H"`, `"RX"`).

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

### Recently done

- ✅ **More builder coverage** (was #1) — `y`, `s`, `t` single-qubit builders
  plus `cz` (controlled-Z) and `cu(gate, control, target)` (arbitrary
  controlled-U) on `Circuit`. Tested in `tests/single_qubit_builders.rs` and
  `tests/controlled_builders.rs`.
- ✅ **Arbitrary controlled-U builder** (was #2) — `mcu(gate, controls, target)`
  and `mcx(controls, target)` on `Circuit`; `toffoli` now delegates to `mcx`.
  Tested in `tests/controlled_builders.rs` (incl. a 3-control C3X).
- ✅ **Circuit diagram printing** (was #1) — `Circuit::diagram()` and a
  `Display` impl render an ASCII diagram (one column per op, `*` controls,
  `|` connectors, `x` for SWAP). Needed a `label` on each `Op`. Tested in
  `tests/diagram.rs`; the demo in `main.rs` prints the circuit.
- ✅ **Richer CLI** (was #1) — a text program format (`src/program.rs`:
  `qubits N`, gate lines, optional `sample SHOTS SEED`, `#` comments, pi-form
  angles) parsed into a `Program`. `main.rs` runs a file, stdin (`-`), or the
  built-in demo, with `--help`. Tested in `tests/program.rs`; sample program
  at `programs/ghz.qsim`. This also covers a GHZ end-to-end path.

- ✅ **GHZ + multi-qubit fixtures** — `tests/ghz.rs` asserts 3- and 4-qubit
  GHZ probabilities and sampling (only 000…/111… outcomes). Salvaged from the
  superseded PR #2 branch, which was an older reimplementation of builders/
  Display that `main` already had.
- ✅ **OpenQASM 2.0 import** — `src/qasm.rs` parses a hand-written subset
  (header/include, one or more `qreg`, the core gate set incl. `cx/cz/ccx/
  swap` and `rx/ry/rz`, `//` and `/* */` comments; `creg/barrier/measure`
  ignored; unsupported features error out). Wired into the CLI by `.qasm`
  extension / `OPENQASM` header. Oracle-tested in `tests/qasm.rs` against the
  equivalent builder circuits; sample at `programs/bell.qasm`.

- ✅ **OpenQASM 2.0 export** — `Circuit::to_qasm()` emits a circuit back to
  OpenQASM (header + `qreg` + one gate per line), lossless for the supported
  subset and round-tripping through `qasm::parse` (rotation angles at full
  `f64` precision). Gates outside the subset (arbitrary controlled-U, C³X)
  return an export error. To support it, `Op::Single` now records `param`
  (rotation angle). CLI: `--emit-qasm`. Tested in `tests/qasm_export.rs`.

- ✅ **More gates (part 1)** — `sdg` (S†), `tdg` (T†), and the phase gate
  `p(λ)` (OpenQASM `u1`). Added to `gates.rs` + builders + both the native
  program format and the QASM import/export maps. Tested in
  `tests/single_qubit_builders.rs`, `tests/program.rs`, `tests/qasm*.rs`, and
  `gates.rs` unit tests.
- ✅ **More gates (part 2a: general single-qubit)** — `u3(θ,φ,λ)` and
  `u2(φ,λ)`. Reshaped `Op::Single`'s param slot from `Option<f64>` to
  `params: Vec<f64>` (0–3 angles); export joins them via `format_params`. Full
  import/export/native-format support, oracle-tested.
- ✅ **More gates (part 2b: controlled rotations)** — `crz(θ)` and controlled
  phase `cp(λ)` (OpenQASM `cu1`). Added `params` to `Op::Controlled` mirroring
  the `Op::Single` reshape; export emits `crz`/`cu1`. Full import/export/
  native-format support, oracle-tested. The QASM subset is now `u2/u3/u1`,
  `rx/ry/rz`, `sdg/tdg`, `cx/cz/crz/cu1/swap/ccx`.

### Suggested next steps (v0.3, not yet started)

Roughly in priority order; none of these exist yet:

1. **Benchmarks** — a committed harness timing gate application across qubit
   counts (state a fixed protocol; no perf claims without it). Only meaningful
   alongside real optimization work (blocked/branch-free kernel, OpenMP-style
   parallelism, or a sparse fast path).
2. **`cu3` / arbitrary controlled-U in QASM** — export currently errors on
   `cu`/`mcu`; a `cu3(θ,φ,λ)` emit + import would round-trip them.
3. **Registers by name in export** — export currently uses a single flat
   `qreg q[n]`; preserving original register names would be nicer (cosmetic).

When adding a gate: add the matrix in `gates.rs`, a builder in `circuit.rs`
(+ `Op` variant and `run` arm if it needs new state machinery), and both a
unit test (matrix/truth-table) and an integration test (circuit behavior).
