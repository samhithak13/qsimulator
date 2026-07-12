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

Both are done in place, so memory is a single `2^n` vector.

## Measurement (planned for v0.1 completion)

Sampling a computational-basis outcome uses the Born rule:
`p(i) = |amplitude(i)|^2`. Post-measurement collapse renormalizes the
surviving amplitudes. Deterministic seeding will be added for
reproducible tests.

## Testing strategy

Every gate has a known truth table / matrix; integration tests assert
exact probabilities for canonical circuits (Bell state, GHZ, single-gate
flips) within a tight epsilon.
