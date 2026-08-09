#!/usr/bin/env python3
"""Cross-validate qsimulator against Qiskit over the OpenQASM 2.0 bridge.

Two phases, each run for `--trials` random circuits:

  * **gates** — generate a random OpenQASM program from the gate set both
    engines implement, and run it through qsimulator (via its `--statevector`
    CLI, which prints the final amplitudes as JSON) and Qiskit's reference
    `Statevector`. This checks that the two engines agree gate for gate.

  * **export** — generate a random program in qsimulator's native text format,
    including multi-controlled gates that OpenQASM 2 has no way to write
    directly. qsimulator runs it, while Qiskit runs the program as qsimulator
    exports it (`--emit-qasm`), i.e. after decomposition. This checks that the
    decomposition an export goes through means the same thing to another tool.

State vectors are compared up to global phase (fidelity), which is the
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
ONE_QUBIT = ["id", "h", "x", "y", "z", "s", "t", "sdg", "tdg"]
ONE_QUBIT_1ANGLE = ["rx", "ry", "rz", "u1"]
TWO_QUBIT = ["cx", "cy", "cz", "ch", "swap"]
TWO_QUBIT_1ANGLE = ["crz", "cu1"]
THREE_QUBIT = ["ccx", "cswap"]


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
            pick = rng.random()
            if pick < 0.45:
                g = rng.choice(TWO_QUBIT)
                lines.append(f"{g} q[{a}],q[{b}];")
            elif pick < 0.8:
                g = rng.choice(TWO_QUBIT_1ANGLE)
                lines.append(f"{g}({angle(rng)}) q[{a}],q[{b}];")
            else:
                lines.append(
                    f"cu3({angle(rng)},{angle(rng)},{angle(rng)}) q[{a}],q[{b}];"
                )
        elif n_qubits >= 3:
            a, b, c = rng.sample(range(n_qubits), 3)
            g = rng.choice(THREE_QUBIT)
            lines.append(f"{g} q[{a}],q[{b}],q[{c}];")
        else:
            q = rng.randrange(n_qubits)
            lines.append(f"x q[{q}];")
    return "\n".join(lines) + "\n"


def random_program(n_qubits: int, n_gates: int, rng: random.Random) -> str:
    """A random circuit in qsimulator's native text format.

    Weighted towards the multi-controlled instructions, since these are the
    ones whose OpenQASM form is a decomposition rather than a single gate; the
    rest are there to hand them a thoroughly entangled input state.
    """
    lines = [f"qubits {n_qubits}"]
    for _ in range(n_gates):
        kind = rng.random()
        if kind < 0.3:
            g = rng.choice(ONE_QUBIT)
            lines.append(f"{g} {rng.randrange(n_qubits)}")
        elif kind < 0.45:
            lines.append(
                f"u3 {angle(rng)} {angle(rng)} {angle(rng)} {rng.randrange(n_qubits)}"
            )
        elif n_qubits >= 2 and kind < 0.55:
            a, b = rng.sample(range(n_qubits), 2)
            lines.append(f"{rng.choice(['cnot', 'cy', 'cz', 'ch', 'swap'])} {a} {b}")
        elif n_qubits >= 2:
            # A multi-controlled gate over a random subset of the register.
            # One control shy of the full register exercises the borrowed-qubit
            # decomposition; the full register exercises the square-root one.
            width = rng.randint(2, n_qubits)
            qubits = rng.sample(range(n_qubits), width)
            operands = " ".join(str(q) for q in qubits)
            if rng.random() < 0.5:
                lines.append(f"mcx {operands}")
            else:
                # theta = 0 is a diagonal gate, which exports by a separate
                # (phase-only) path, so hit it deliberately every so often.
                theta = "0" if rng.random() < 0.3 else angle(rng)
                lines.append(f"mcu3 {theta} {angle(rng)} {angle(rng)} {operands}")
        else:
            lines.append(f"x {rng.randrange(n_qubits)}")
    return "\n".join(lines) + "\n"


def qsim(args: list[str], program: str, binary: str) -> str:
    result = subprocess.run(
        [binary, *args, "-"],
        input=program,
        capture_output=True,
        text=True,
        check=True,
    )
    return result.stdout


def qsim_statevector(qasm: str, binary: str) -> np.ndarray:
    data = json.loads(qsim(["--statevector"], qasm, binary))
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


def compare(
    label: str,
    trial: int,
    source: str,
    a: np.ndarray,
    b: np.ndarray,
    tol: float,
) -> float | None:
    """Fidelity of the two state vectors, or None (after printing) on a
    mismatch."""
    if a.shape != b.shape:
        print(f"FAIL ({label} trial {trial}): shape {a.shape} vs {b.shape}\n{source}")
        return None
    f = fidelity(a, b)
    if 1.0 - f > tol:
        print(f"FAIL ({label} trial {trial}): fidelity {f:.3e} below 1 - {tol:g}")
        print(source)
        return None
    return f


def run_gate_phase(args, binary: str, rng: random.Random) -> float | None:
    """Both engines run the same OpenQASM program."""
    worst = 1.0
    for trial in range(args.trials):
        n_qubits = rng.randint(1, args.max_qubits)
        n_gates = rng.randint(2 * n_qubits, 6 * n_qubits)
        qasm = random_qasm(n_qubits, n_gates, rng)

        f = compare(
            "gates",
            trial,
            qasm,
            qsim_statevector(qasm, binary),
            qiskit_statevector(qasm),
            args.tol,
        )
        if f is None:
            return None
        worst = min(worst, f)
    return worst


def run_export_phase(args, binary: str, rng: random.Random) -> float | None:
    """qsimulator runs a native program; Qiskit runs what qsimulator exports,
    so any error in the export decomposition shows up as a mismatch."""
    worst = 1.0
    for trial in range(args.trials):
        n_qubits = rng.randint(2, args.max_qubits)
        n_gates = rng.randint(n_qubits, 3 * n_qubits)
        program = random_program(n_qubits, n_gates, rng)
        exported = qsim(["--emit-qasm"], program, binary)

        f = compare(
            "export",
            trial,
            f"{program}\n--- exported as ---\n{exported}",
            qsim_statevector(program, binary),
            qiskit_statevector(exported),
            args.tol,
        )
        if f is None:
            return None
        worst = min(worst, f)
    return worst


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trials", type=int, default=500, help="trials per phase")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--tol", type=float, default=1e-9)
    parser.add_argument("--max-qubits", type=int, default=5)
    parser.add_argument("--binary", default=None, help="path to the qsimulator binary")
    args = parser.parse_args()

    binary = args.binary or build_binary()
    rng = random.Random(args.seed)

    for label, phase in (("gates", run_gate_phase), ("export", run_export_phase)):
        worst = phase(args, binary, rng)
        if worst is None:
            return 1
        print(
            f"OK: {args.trials} {label} trials agree with Qiskit "
            f"(worst fidelity {worst:.15f}, tol {args.tol:g})"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
