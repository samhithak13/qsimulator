#!/usr/bin/env python3
"""Cross-validate qsimulator against Qiskit over the OpenQASM 2.0 bridge.

For each trial this generates a random circuit from the gate set both engines
implement, emits it as OpenQASM 2.0, and runs it through:

  * qsimulator, via its `--statevector` CLI (final amplitudes as JSON), and
  * Qiskit's reference `Statevector`.

The two state vectors are compared up to global phase (fidelity), which is the
physically meaningful notion of state equality and is robust to per-gate
global-phase conventions. A single mismatch prints the offending program and
exits non-zero.

Usage:
    python crossval/compare.py [--trials N] [--seed S] [--tol T] [--binary PATH]
"""

from __future__ import annotations

import argparse
import json
import os
import random
import subprocess
import sys

import numpy as np

# Gates shared by qsimulator and Qiskit's qelib1, grouped by (qubits, angles).
ONE_QUBIT = ["h", "x", "y", "z", "s", "t", "sdg", "tdg"]
ONE_QUBIT_1ANGLE = ["rx", "ry", "rz", "u1"]
TWO_QUBIT = ["cx", "cz", "swap"]
TWO_QUBIT_1ANGLE = ["crz", "cu1"]


def angle(rng: random.Random) -> str:
    # repr() round-trips exactly through both parsers.
    return repr(rng.uniform(-2.0 * np.pi, 2.0 * np.pi))


def random_qasm(n_qubits: int, n_gates: int, rng: random.Random) -> str:
    lines = ["OPENQASM 2.0;", 'include "qelib1.inc";', f"qreg q[{n_qubits}];"]
    for _ in range(n_gates):
        kind = rng.random()
        if kind < 0.35:
            g = rng.choice(ONE_QUBIT)
            q = rng.randrange(n_qubits)
            lines.append(f"{g} q[{q}];")
        elif kind < 0.55:
            g = rng.choice(ONE_QUBIT_1ANGLE)
            q = rng.randrange(n_qubits)
            lines.append(f"{g}({angle(rng)}) q[{q}];")
        elif kind < 0.70:
            q = rng.randrange(n_qubits)
            lines.append(f"u2({angle(rng)},{angle(rng)}) q[{q}];")
        elif kind < 0.82:
            q = rng.randrange(n_qubits)
            lines.append(f"u3({angle(rng)},{angle(rng)},{angle(rng)}) q[{q}];")
        elif n_qubits >= 2 and kind < 0.94:
            a, b = rng.sample(range(n_qubits), 2)
            if rng.random() < 0.5:
                g = rng.choice(TWO_QUBIT)
                lines.append(f"{g} q[{a}],q[{b}];")
            else:
                g = rng.choice(TWO_QUBIT_1ANGLE)
                lines.append(f"{g}({angle(rng)}) q[{a}],q[{b}];")
        elif n_qubits >= 3:
            a, b, c = rng.sample(range(n_qubits), 3)
            lines.append(f"ccx q[{a}],q[{b}],q[{c}];")
        else:
            q = rng.randrange(n_qubits)
            lines.append(f"x q[{q}];")
    return "\n".join(lines) + "\n"


def qsim_statevector(qasm: str, binary: str) -> np.ndarray:
    result = subprocess.run(
        [binary, "--statevector", "-"],
        input=qasm,
        capture_output=True,
        text=True,
        check=True,
    )
    data = json.loads(result.stdout)
    return np.array([complex(re, im) for re, im in data], dtype=complex)


def qiskit_statevector(qasm: str) -> np.ndarray:
    from qiskit import qasm2
    from qiskit.quantum_info import Statevector

    circuit = qasm2.loads(qasm, custom_instructions=qasm2.LEGACY_CUSTOM_INSTRUCTIONS)
    # Qiskit and qsimulator both index little-endian (qubit 0 is the least
    # significant bit), so no reordering is needed.
    return np.asarray(Statevector(circuit).data, dtype=complex)


def fidelity(a: np.ndarray, b: np.ndarray) -> float:
    """State fidelity |<a|b>|^2, which ignores an overall global phase."""
    return float(np.abs(np.vdot(a, b)) ** 2)


def build_binary() -> str:
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    subprocess.run(["cargo", "build", "--quiet"], cwd=root, check=True)
    return os.path.join(root, "target", "debug", "qsimulator")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trials", type=int, default=500)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--tol", type=float, default=1e-9)
    parser.add_argument("--max-qubits", type=int, default=5)
    parser.add_argument("--binary", default=None, help="path to the qsimulator binary")
    args = parser.parse_args()

    binary = args.binary or build_binary()
    rng = random.Random(args.seed)

    worst = 1.0
    for trial in range(args.trials):
        n_qubits = rng.randint(1, args.max_qubits)
        n_gates = rng.randint(2 * n_qubits, 6 * n_qubits)
        qasm = random_qasm(n_qubits, n_gates, rng)

        a = qsim_statevector(qasm, binary)
        b = qiskit_statevector(qasm)

        if a.shape != b.shape:
            print(f"FAIL (trial {trial}): shape {a.shape} vs {b.shape}\n{qasm}")
            return 1

        f = fidelity(a, b)
        worst = min(worst, f)
        if 1.0 - f > args.tol:
            print(f"FAIL (trial {trial}): fidelity {f:.3e} below 1 - {args.tol:g}")
            print(qasm)
            return 1

    print(
        f"OK: {args.trials} trials agree with Qiskit "
        f"(worst fidelity {worst:.15f}, tol {args.tol:g})"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
