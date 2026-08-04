// Bell state in OpenQASM 2.0.
// Run with:  cargo run -- programs/bell.qasm
OPENQASM 2.0;
include "qelib1.inc";
qreg q[2];
creg c[2];
h q[0];
cx q[0],q[1];
measure q -> c;
