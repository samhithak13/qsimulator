# Cross-validation against Qiskit

`compare.py` checks qsimulator against [Qiskit](https://www.ibm.com/quantum/qiskit)
as an independent reference implementation. Each trial builds a random circuit
from the gate set both engines share, emits it as OpenQASM 2.0, runs it through
qsimulator (`--statevector`) and Qiskit's `Statevector`, and compares the two
state vectors up to global phase.

## Running

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r crossval/requirements.txt
python crossval/compare.py --trials 500
```

The script builds the qsimulator binary itself and drives it through
`qsimulator --statevector -`. Options: `--trials`, `--seed`, `--tol`,
`--max-qubits`, and `--binary` (to point at a prebuilt binary).

## What it establishes

Over randomized circuits spanning the shared gate set — H, X, Y, Z, S, T and
their daggers, `u1`/`u2`/`u3`, `rx`/`ry`/`rz`, `cx`, `cz`, `crz`, `cu1`,
`swap`, `ccx` — qsimulator's amplitudes match Qiskit's to floating-point
precision (state fidelity within `1e-9` of 1).

Comparison is done up to a global phase, via `|<a|b>|^2`. That is the
physically meaningful notion of state equality and avoids spurious failures
from per-gate global-phase conventions, which differ between implementations.
Both engines index basis states little-endian (qubit 0 is the least
significant bit), so no qubit reordering is applied.
