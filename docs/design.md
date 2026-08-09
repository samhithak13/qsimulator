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

Both are done in place, so memory is a single `2^n` vector. `apply_1q` walks
the amplitude vector as `chunks_exact_mut(2·2^target)` and splits each block
into the target bit's `|0>`/`|1>` halves, so the inner loop is bounds-check
free; the gate-matrix entries are hoisted into locals in every kernel.

The controlled and multi-controlled kernels share one private `apply_masked`
routine: the same target-bit walk with a per-pair test that the control bits
are set (a single control being just a one-bit mask).

With the optional `parallel` feature, `apply_1q` and `apply_masked` run across
threads via rayon: they parallelize over blocks for a low target qubit (coarse
tiles, to keep tasks large) and within a block's two halves for a high one.
The feature is off by default, so the default build keeps its single
`num-complex` dependency and the `forbid(unsafe_code)` guarantee (rayon's safe
parallel iterators need no unsafe here). See `benches/` for the harness and
numbers.

The gate set covers the Paulis (X, Y, Z), Hadamard, the phase gates S and
T (and their daggers), the phase gate `p(λ)`, the general single-qubit
`u2`/`u3`, and the continuous rotations `rx(θ)`, `ry(θ)`, `rz(θ)`. Each
rotation is the standard `exp(-i·θ/2·P)` for its Pauli `P`, so e.g. `rx(π)`
equals X up to the global phase `-i`.

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

Every quantity is checked against an independent oracle rather than a
hand-copied expected value. Gate matrices are asserted against their
algebraic identities; canonical circuits (Bell, GHZ, single-gate flips)
are asserted to exact probabilities within a tight epsilon; the OpenQASM
importer is validated against the equivalent builder circuit and the
exporter by round-tripping back through it.

Beyond the in-tree tests, `crossval/compare.py` validates the engine against
Qiskit: it runs random circuits over the shared gate set through both
implementations (via OpenQASM) and compares the state vectors up to global
phase. The `--statevector` CLI flag exposes the final amplitudes as JSON for
this harness.

## Status

Implemented and covered by tests:

- **State vector** (`State`) — construction, amplitudes, per-index
  probability, norm; single-, controlled-, and multi-controlled gate
  application; qubit swap.
- **Gates** (`gates`) — the identity, X, Y, Z, H, S, T, S†, T†, phase `p(λ)`,
  `u2`, `u3`, and `rx`/`ry`/`rz`.
- **Measurement** — `prob_qubit_one`, `measure_qubit`, `measure_all`, and
  `Circuit::sample`, backed by the seedable `rng::Rng`.
- **Circuit** (`circuit::Circuit`) — builder methods for the full gate set,
  including `cnot`, `cz`, `crz`, `cp`, `cu3`, `swap`, `toffoli`, and the
  general `cu`/`mcx`/`mcu`; execution; ASCII diagrams; OpenQASM export at any
  control width.
- **Front ends** — a text program parser (`program`), an OpenQASM 2.0
  importer (`qasm`), and a CLI (`main`).

Over 150 unit and integration tests at ~98% line coverage; CI runs fmt,
clippy (`-D warnings`), a warning-clean `cargo doc`, build, test, the Qiskit
cross-validation, parser fuzzing, `cargo audit`, and a coverage floor.

## Module layout

- `src/state.rs` — the state vector, gate application (`apply_1q`,
  `apply_controlled_1q`, `apply_multi_controlled_1q`), measurement, and
  `swap_qubits`.
- `src/gates.rs` — the 2×2 gate matrices (`type Gate = [[Complex64; 2]; 2]`)
  and their unit tests.
- `src/circuit.rs` — the `Op` enum, the `Circuit` builder, `run`, `sample`,
  `diagram`, and `to_qasm`.
- `src/program.rs` — the text program parser (`parse`, `Program`,
  `SampleSpec`); `parse_angle` is shared with the QASM importer.
- `src/qasm.rs` — the OpenQASM 2.0 subset importer.
- `src/rng.rs` — the seedable RNG.
- `src/main.rs` — the CLI; dispatches to the QASM importer by `.qasm`
  extension or `OPENQASM` header, otherwise the text format.
- `tests/` — one file per area: measurement, rotations, swap, toffoli, the
  builder groups, diagrams, GHZ, and QASM import/export.

`Op` variants carry a `label` and `params` used by `diagram()` and
`to_qasm()`; neither affects execution. A new gate sets both.

## Conventions

These are load-bearing; changing them breaks results silently.

- **Little-endian qubit order** — bit `q` of a basis index is qubit `q`. The
  2-qubit index `0b10` means qubit 1 is |1> and qubit 0 is |0>.
- **Gate storage** is row-major `[[m00, m01], [m10, m11]]` over `Complex64`.
- **Rotations** are `exp(-i·θ/2·P)`, so `R_axis(π)` equals its Pauli up to the
  global phase `-i` (the tests assert this).
- **RNG output is a pure function of the seed.** Reproducibility depends on
  the xorshift64 recurrence and SplitMix64 seeding staying fixed.
- Controls must differ from their target, and each amplitude pair is touched
  once, from the target-bit-0 side.

## Development

```bash
cargo fmt --all -- --check
cargo clippy --all-targets -- -D warnings
cargo build
cargo test
```

New functionality lands with its test in the same commit, and all four checks
pass before committing. Adding a gate touches its matrix in `gates.rs`, a
builder in `circuit.rs` (and an `Op` arm if it needs new state machinery), the
two front-end maps if it should be importable, and both a matrix/truth-table
unit test and a circuit-level integration test.

An arbitrary controlled-U (`cu`) has no direct OpenQASM 2 gate, so `to_qasm`
decomposes it into a control phase (`u1`) plus `cu3` via a single-qubit Euler
decomposition (`gates::u3_decompose`, which writes any 2×2 unitary as
`e^{iγ}·u3(θ,φ,λ)`). The `cu3` convention is cross-validated against Qiskit.

### Multi-controlled export

OpenQASM 2 stops at `ccx`, so wider multi-controlled gates are decomposed
(`emit_mcx`, `emit_mcu`, `emit_mcphase` in `circuit.rs`). No ancillas are ever
allocated — the register is exactly what the user declared.

`C^m(X)` uses one of two Barenco constructions, chosen by whether the register
has a qubit the operation does not touch:

- **A spare qubit exists** — it is *borrowed*: its state is unknown, so the
  controls are split in half, the borrowed qubit is toggled by the first half
  and used as an extra control for the second, and the pair is run twice. The
  unknown initial value cancels out of the target and the borrowed qubit is
  left exactly as it was found (the gates are permutations, so this holds on
  superpositions too). Each half recurses, giving `O(m²)` Toffolis.
- **The controls plus the target are the whole register** — the square-root
  recursion, with `V·V = X`: `C^m(X) = C(V)·C^{m-1}(X)·C(V†)·C^{m-1}(X)·
  C^{m-1}(V)`. Every inner operation now leaves a qubit untouched, so it lands
  in the borrowed-qubit case.

`C^m(U)` writes `U` as `e^{iγ''}·A·X·B·X·C` with `A·B·C = I` (`A`, `B`, `C`
being `rz`/`ry` products of the Euler angles). Only the two `X`s and the phase
`γ''` are conditioned on the controls: off the all-controls-set subspace,
`A·B·C` collapses to the identity. A diagonal `U` (multi-controlled Z, S, T,
phase) skips the `X`s and exports as two phase terms. Multi-controlled phases
recurse through the same identity that decomposes `cu1`.

The decomposition is exact, not up-to-global-phase: `tests/qasm_export.rs`
compares re-imported amplitudes elementwise. The single exception is an `mcu`
with *no* controls, where the phase is genuinely global and OpenQASM 2 has no
way to write it.

## Roadmap

- A SIMD or blocked-complex-arithmetic fast path, if the memory-bound ceiling
  ever needs raising.
