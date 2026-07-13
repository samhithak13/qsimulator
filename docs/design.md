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
