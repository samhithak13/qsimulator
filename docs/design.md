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
- `Circuit::sample(shots, seed)` — a histogram over `shots` measurements. When
  nothing in the circuit branches it runs once and measures clones of the final
  state; when a collapse can affect what follows — a mid-circuit measurement, a
  reset, or a noise channel — it re-runs the circuit per shot, since there is
  no single final state to draw from.

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
Qiskit in four phases. Two compare state vectors up to global phase — random
circuits over the shared gate set, and qsimulator's own *exported* form read
back by Qiskit, which is what checks the decompositions. The other two compare
sampled distributions against Qiskit Aer, because a circuit that collapses or
decoheres has no single state vector: one covers mid-circuit measurement,
reset and classical feed-forward, the other noise channels.

Each sampled phase was checked against the bug it exists to catch, not just
for passing: the measurement phase originally agreed with Aer even with the
collapse ignored entirely, because a collapse commutes with anything diagonal
or permutation-like in the computational basis. Its generator now sandwiches
the collapse between basis-changing gates. See `crossval/README.md`.

## Status

Implemented and covered by tests:

- **State vector** (`State`) — construction, amplitudes, per-index
  probability, norm; single-, controlled-, and multi-controlled gate
  application; qubit swap.
- **Gates** (`gates`) — the identity, X, Y, Z, H, S, T, S†, T†, `sx`/`sxdg`,
  phase `p(λ)`, `u2`, `u3`, and `rx`/`ry`/`rz`.
- **Measurement** — `prob_qubit_one`, `measure_qubit`, `measure_all`, and
  `Circuit::sample`, backed by the seedable `rng::Rng`.
- **Non-unitary operations** — mid-circuit `measure` and `reset`, classical
  feed-forward (`if_classical_eq`) over one classical register, and noise
  channels via `apply_kraus`. Each makes execution depend on the RNG stream,
  which `run_seeded` chooses.
- **Noise** (`noise`) — depolarizing, bit flip, phase flip, amplitude damping
  and phase damping, with a trace-preservation check.
- **Circuit** (`circuit::Circuit`) — builder methods for the full gate set,
  including `cnot`, `cz`, `crz`, `cp`, `cu3`, `swap`, `toffoli`, and the
  general `cu`/`mcx`/`mcu`; execution; ASCII diagrams; OpenQASM export at any
  control width.
- **Front ends** — a text program parser (`program`), an OpenQASM 2.0
  importer (`qasm`) covering `gate` declarations and all of `qelib1`, and a
  CLI (`main`).

Over 190 unit and integration tests at ~97% line coverage; CI runs fmt,
clippy (`-D warnings`), a warning-clean `cargo doc`, build, test, the Qiskit
cross-validation, parser fuzzing, `cargo audit`, and a coverage floor.

## Module layout

- `src/state.rs` — the state vector, gate application (`apply_1q`,
  `apply_controlled_1q`, `apply_multi_controlled_1q`), measurement and
  collapse, `apply_kraus` for noise, and `swap_qubits`.
- `src/gates.rs` — the 2×2 gate matrices (`type Gate = [[Complex64; 2]; 2]`)
  and their unit tests.
- `src/circuit.rs` — the `Op` enum, the `Circuit` builder, `run`, `sample`,
  `diagram`, and `to_qasm`.
- `src/program.rs` — the text program parser (`parse`, `Program`,
  `SampleSpec`); `parse_angle` is shared with the QASM importer.
- `src/qasm.rs` — the OpenQASM 2.0 subset importer.
- `src/expr.rs` — the angle expression evaluator, shared by both front ends;
  `gate` bodies need arithmetic over their formal parameters.
- `src/noise.rs` — the standard single-qubit channels as Kraus operators, and
  the trace-preservation check that rejects a non-physical one.
- `src/rng.rs` — the seedable RNG.
- `src/main.rs` — the CLI; dispatches to the QASM importer by `.qasm`
  extension or `OPENQASM` header, otherwise the text format.
- `tests/` — one file per area: measurement and collapse, noise, rotations,
  swap, toffoli, the builder groups, diagrams, GHZ, QASM import/export, the
  text program format, the CLI, and property and robustness tests.

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

## Measurement and classical control

`Circuit` is no longer a sequence of unitaries. Three operations collapse the
state, and all three make execution depend on the RNG stream, which
`run_seeded` chooses:

- `measure(q)` collapses `q` and records the outcome in the circuit's one
  classical register, which is as wide as the quantum register — bit `i` is the
  last outcome written to bit `i`, and unmeasured bits read 0. `measure_into`
  picks a different destination bit, which an imported program needs when it
  compacts measurements into low bits.
- `reset(q)` collapses `q` and forces it to |0>.
- `if_classical_eq(value, ...)` runs a block of gates only when the whole
  classical register equals `value` — OpenQASM's `if`.

Two consequences are load-bearing:

- **Trailing measurements are readout.** Nothing follows them, so collapsing
  would discard the prepared state and report one arbitrary branch, and since
  nearly every written-out program ends in a measurement that would make `run`
  seed-dependent for almost every real file. Sampling draws the same
  distribution either way. A trailing *reset* is not readout: it changes the
  state that gets sampled.
- **`sample` re-runs the circuit per shot** when a collapse can affect what
  follows, because there is no single final state to draw from. Otherwise it
  keeps the cheap path of running once and measuring clones.

A conditional block holds only gates. A measurement inside one would change the
value being tested part-way through, so guarding each statement — all
OpenQASM's single-statement `if` can express — would stop matching the block as
a whole; a nested conditional has no OpenQASM form at all. Both panic rather
than being dropped.

None of this can be checked by comparing state vectors, since a branching
circuit has none. The third `crossval` phase samples against Qiskit Aer
instead; see `crossval/README.md` for why its generator is shaped the way it is.

## Noise

A noise channel maps a density matrix to `sum_i K_i rho K_i^dagger`. This
engine stores a state vector, so channels are simulated by the
**quantum-trajectory** method instead: each shot samples one Kraus operator
`K_i` with probability `<psi|K_i^dagger K_i|psi>` and renormalizes
(`State::apply_kraus`). Averaged over shots the ensemble reproduces the density
matrix.

This fits the existing machinery rather than needing new machinery: `sample`
already re-runs the circuit per shot with one RNG stream, which is exactly the
structure trajectories need. The cost is statistical — a single `run` is one
trajectory, not the average, and results converge as `1/sqrt(shots)`. A density
matrix would give exact answers in one pass at the price of squaring the memory,
from `2^n` amplitudes to `2^n x 2^n`.

Branch probabilities are computed as `<psi|K^dagger K|psi>` in a single
read-only pass (`State::expectation_1q`), so a channel costs no clone of the
state vector — which the obvious implementation, applying each `K_i` to a copy
to measure its norm, would need once per operator.

Channels are rejected at build time unless they are trace preserving
(`noise::is_trace_preserving`), since a set that is not would silently change
the total probability. Noise also makes a circuit stochastic in the same way a
reset does: it always branches, even as the final operation.

OpenQASM 2 has no syntax for a channel, so exporting a noisy circuit returns
`ExportError::Noise` rather than writing out a circuit that quietly differs
from the one that ran. That is also why the noise cross-validation phase builds
its Qiskit circuit directly instead of sharing a source file.

## Density matrices

`DensityMatrix` is the second backend. Where a state vector must sample one
Kraus operator per shot, `rho` carries the mixture itself, so a channel is
applied exactly: `rho -> sum_i K_i rho K_i^dagger` in one step, with no
sampling error and no shots. `purity` (`Tr(rho^2)`) then measures directly how
much the noise cost — 1 for a pure state, `1/2^n` for the maximally mixed one.

Unitary evolution is conjugation, `rho -> U rho U^dagger`, done as two sweeps
over the same butterfly: `U` along the row index, then `conj(U)` along the
column index. That reuses one routine instead of materializing `U rho` and
multiplying again.

The cost is `4^n` entries against `2^n`, so the ceiling is
`MAX_DENSITY_QUBITS` = 12 (~268 MB) rather than 30. Twelve qubits here is the
same memory as twenty-four there.

Two boundaries are worth stating:

- **An unread measurement is exact.** It is precisely the channel
  `rho -> P_0 rho P_0 + P_1 rho P_1`, which erases coherence between outcomes
  and is deterministic. So is `reset`. Neither needs sampling here.
- **Classical feed-forward is not representable.** Branching on an outcome
  needs a distribution over classical registers, each with its own matrix;
  `rho` alone has no classical state. `run_density` returns
  `DensityError::ClassicalFeedForward` rather than approximating it, and such
  a circuit should be sampled instead.

Because both backends are exact for unitaries and both implement the same
channels, they check each other: the tests assert that trajectory sampling
converges to what the density matrix says exactly.

## Roadmap

Nothing outstanding. The OpenQASM bridge is complete in both directions: every
builder gate exports, and every gate Qiskit emits — including `gate`
declarations and the whole of `qelib1` — imports. `measure`, `reset` and `if`
are all honoured, leaving only `opaque`, which declares no body to simulate.

Both noise backends now exist — trajectories for reach, density matrices for
exactness — so the remaining direction is OpenQASM 3, or making feed-forward
representable in the density backend by carrying a classical mixture.

### Not planned: a SIMD fast path

Earlier notes listed a SIMD or blocked-complex-arithmetic kernel "if the
memory-bound ceiling ever needs raising". It does not, and the reason is
arithmetic intensity rather than effort.

Applying a 2x2 gate to one amplitude pair moves 64 bytes (two `Complex64`
loaded, two stored) to do four complex multiplies and two complex adds — about
28 flops, so roughly 0.4 flops per byte. That is one to two orders of magnitude
below the flops-per-byte a core can sustain against DRAM, so the kernel waits
on memory, not on the multiplier. Widening the arithmetic makes the part that
was never the bottleneck faster.

The `parallel` benchmarks corroborate it: eight cores buy 1.2x–2.4x, not 8x,
and the best case is exactly the one that streams two large contiguous halves.
That shape is bandwidth-limited scaling.

There is also a cost. Portable SIMD (`std::simd`) is nightly-only, and stable
SIMD means `core::arch` intrinsics, which are `unsafe` — so the fast path would
trade away `#![forbid(unsafe_code)]`, a guarantee worth more than a few percent
on a memory-bound loop. The kernels are already shaped for the compiler to
auto-vectorize what it can: a bounds-check-free walk over `chunks_exact_mut`
with the gate entries hoisted into locals.

Should this ever be revisited, the profitable direction is reducing traffic
rather than widening arithmetic — fusing adjacent single-qubit gates on the
same target so one sweep does the work of several, or blocking several gates
over a cache-resident slice of the vector.
