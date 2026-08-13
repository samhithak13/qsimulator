//! Integration tests for the OpenQASM 2.0 subset importer.
//!
//! The oracle is the builder API: an imported circuit must produce the same
//! state as the equivalent hand-built circuit.

use approx::assert_relative_eq;
use qsimulator::{qasm, Circuit, State};

/// Assert two states agree amplitude by amplitude — stricter than comparing
/// probabilities, since an expanded `gate` body must reproduce phases too.
fn assert_same_amplitudes(a: &State, b: &State) {
    assert_eq!(a.n_qubits(), b.n_qubits());
    for (i, (x, y)) in a.amplitudes().iter().zip(b.amplitudes()).enumerate() {
        assert!((x - y).norm() < 1e-12, "amplitude {i} differs: {x} vs {y}");
    }
}

/// Assert two states have the same per-basis-state probabilities.
fn assert_same_probs(a: &State, b: &State) {
    assert_eq!(a.n_qubits(), b.n_qubits());
    for i in 0..a.amplitudes().len() {
        assert_relative_eq!(a.probability(i), b.probability(i), epsilon = 1e-12);
    }
}

#[test]
fn imports_bell_state() {
    let src = "\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[2];
creg c[2];
h q[0];
cx q[0],q[1];
measure q -> c;
";
    let imported = qasm::parse(src).expect("should parse").run();

    let mut expected = Circuit::new(2);
    expected.h(0).cnot(0, 1);
    assert_same_probs(&imported, &expected.run());
}

#[test]
fn imports_ghz_and_ignores_barrier() {
    let src = "\
OPENQASM 2.0;
qreg q[3];
h q[0];
barrier q;
cx q[0],q[1];
cx q[1],q[2];
";
    let state = qasm::parse(src).expect("should parse").run();
    assert_relative_eq!(state.probability(0b000), 0.5, epsilon = 1e-12);
    assert_relative_eq!(state.probability(0b111), 0.5, epsilon = 1e-12);
}

#[test]
fn imports_rotation_with_pi_angle() {
    // rx(pi) flips |0> to |1> up to global phase.
    let state = qasm::parse("qreg q[1];\nrx(pi) q[0];\n")
        .expect("should parse")
        .run();
    assert_relative_eq!(state.probability(1), 1.0, epsilon = 1e-12);
}

#[test]
fn imports_ccx_toffoli() {
    // Both controls set -> target flips.
    let src = "qreg q[3];\nx q[0];\nx q[1];\nccx q[0],q[1],q[2];\n";
    let state = qasm::parse(src).expect("should parse").run();
    assert_relative_eq!(state.probability(0b111), 1.0, epsilon = 1e-12);
}

#[test]
fn imports_sdg_tdg_and_u1_phase() {
    // x then t then tdg then sdg then s -> net phase +1 on |1> (identity).
    let src = "qreg q[1];\nx q[0];\nt q[0];\ntdg q[0];\nsdg q[0];\ns q[0];\n";
    let imported = qasm::parse(src).expect("should parse").run();
    assert_relative_eq!(imported.amplitudes()[1].re, 1.0, epsilon = 1e-12);
    assert_relative_eq!(imported.amplitudes()[1].im, 0.0, epsilon = 1e-12);

    // u1(pi) is the phase gate; on |1> it matches Z.
    let mut expected = Circuit::new(1);
    expected.x(0).p(std::f64::consts::PI, 0);
    let via_u1 = qasm::parse("qreg q[1];\nx q[0];\nu1(pi) q[0];\n").unwrap();
    assert_same_probs(&via_u1.run(), &expected.run());
    assert_relative_eq!(via_u1.run().amplitudes()[1].re, -1.0, epsilon = 1e-12);

    // `p(lambda)` is accepted as an alias for `u1(lambda)`.
    let via_p = qasm::parse("qreg q[1];\nx q[0];\np(pi) q[0];\n").unwrap();
    assert_same_probs(&via_p.run(), &expected.run());
}

#[test]
fn imports_u2_and_u3() {
    // u3(pi,0,pi) = X, so on |0> gives |1>.
    let s1 = qasm::parse("qreg q[1];\nu3(pi,0,pi) q[0];\n")
        .expect("should parse")
        .run();
    assert_relative_eq!(s1.probability(1), 1.0, epsilon = 1e-12);

    // u2(0,pi) = H, so on |0> gives an equal superposition.
    let s2 = qasm::parse("qreg q[1];\nu2(0,pi) q[0];\n")
        .expect("should parse")
        .run();
    assert_relative_eq!(s2.probability(0), 0.5, epsilon = 1e-12);
    assert_relative_eq!(s2.probability(1), 0.5, epsilon = 1e-12);
}

#[test]
fn error_u3_wrong_angle_count() {
    let err = qasm::parse("qreg q[1];\nu3(pi,0) q[0];\n")
        .unwrap_err()
        .to_string();
    assert!(err.contains("3 angle(s)"), "{err}");
}

#[test]
fn imports_controlled_rotations() {
    // cu1(pi) on |11> negates that amplitude; cp is accepted as an alias.
    for src in [
        "qreg q[2];\nx q[0];\nx q[1];\ncu1(pi) q[0],q[1];\n",
        "qreg q[2];\nx q[0];\nx q[1];\ncp(pi) q[0],q[1];\n",
    ] {
        let state = qasm::parse(src).expect("should parse").run();
        assert_relative_eq!(state.probability(0b11), 1.0, epsilon = 1e-12);
        assert_relative_eq!(state.amplitudes()[0b11].re, -1.0, epsilon = 1e-12);
    }

    // crz round-trips through export.
    let c = qasm::parse("qreg q[2];\ncrz(0.7) q[0],q[1];\n").expect("should parse");
    let reimported = qasm::parse(&c.to_qasm().unwrap()).unwrap();
    assert_same_probs(&reimported.run(), &c.run());

    // cu3(pi, 0, pi) is a controlled-X, so with the control set it flips the
    // target: |01> -> |11>.
    let cnot_like = qasm::parse("qreg q[2];\nx q[0];\ncu3(pi,0,pi) q[0],q[1];\n")
        .expect("should parse")
        .run();
    assert_relative_eq!(cnot_like.probability(0b11), 1.0, epsilon = 1e-12);
}

#[test]
fn imports_cy_ch_cswap() {
    // cy: control set, Y|0> = i|1>, so amplitude of |11> is i.
    let cy = qasm::parse("qreg q[2];\nx q[0];\ncy q[0],q[1];\n")
        .expect("should parse")
        .run();
    assert_relative_eq!(cy.probability(0b11), 1.0, epsilon = 1e-12);
    assert_relative_eq!(cy.amplitudes()[0b11].im, 1.0, epsilon = 1e-12);

    // ch: control set, H|0> spreads the target evenly.
    let ch = qasm::parse("qreg q[2];\nx q[0];\nch q[0],q[1];\n")
        .expect("should parse")
        .run();
    assert_relative_eq!(ch.probability(0b01), 0.5, epsilon = 1e-12);
    assert_relative_eq!(ch.probability(0b11), 0.5, epsilon = 1e-12);

    // cswap: control set swaps the other two (|011> -> |101>); control clear
    // leaves them (|010> stays).
    let swapped = qasm::parse("qreg q[3];\nx q[0];\nx q[1];\ncswap q[0],q[1],q[2];\n")
        .expect("should parse")
        .run();
    assert_relative_eq!(swapped.probability(0b101), 1.0, epsilon = 1e-12);

    let unswapped = qasm::parse("qreg q[3];\nx q[1];\ncswap q[0],q[1],q[2];\n")
        .expect("should parse")
        .run();
    assert_relative_eq!(unswapped.probability(0b010), 1.0, epsilon = 1e-12);
}

#[test]
fn imports_cz_and_swap() {
    let src = "qreg q[2];\nx q[0];\nswap q[0],q[1];\ncz q[0],q[1];\n";
    let imported = qasm::parse(src).expect("should parse").run();

    let mut expected = Circuit::new(2);
    expected.x(0).swap(0, 1).cz(0, 1);
    assert_same_probs(&imported, &expected.run());
}

#[test]
fn multiple_registers_map_to_flat_space() {
    // q[0..2] -> 0,1 ; r[0] -> 2. cx q[0],r[0] entangles qubit 0 and qubit 2.
    let src = "qreg q[2];\nqreg r[1];\nh q[0];\ncx q[0],r[0];\n";
    let imported = qasm::parse(src).expect("should parse").run();

    let mut expected = Circuit::new(3);
    expected.h(0).cnot(0, 2);
    assert_same_probs(&imported, &expected.run());
}

#[test]
fn strips_line_and_block_comments() {
    let src = "\
// leading comment
qreg q[1]; /* inline */ x q[0]; // trailing
/* multi
   line */
";
    let state = qasm::parse(src).expect("should parse").run();
    assert_relative_eq!(state.probability(1), 1.0, epsilon = 1e-12);
}

#[test]
fn error_unsupported_gate() {
    // A name that is not a builtin and not declared in the file. (`sx` used to
    // stand in here, before it became a supported gate.)
    let err = qasm::parse("qreg q[1];\nnotagate q[0];\n")
        .unwrap_err()
        .to_string();
    assert!(err.contains("unsupported gate `notagate`"), "{err}");
}

#[test]
fn error_unsupported_feature() {
    // `opaque` declares a gate with no body, so there is nothing to simulate.
    // It is the last unsupported feature; `reset` and `if` each stood here
    // before they were implemented.
    let err = qasm::parse("qreg q[1];\nopaque mystery a;\n")
        .unwrap_err()
        .to_string();
    assert!(
        err.contains("unsupported OpenQASM feature `opaque`"),
        "{err}"
    );
}

#[test]
fn error_no_qreg() {
    let err = qasm::parse("OPENQASM 2.0;\nh q[0];\n")
        .unwrap_err()
        .to_string();
    assert!(
        err.contains("no `qreg`") || err.contains("unknown register"),
        "{err}"
    );
}

#[test]
fn error_qubit_out_of_range() {
    let err = qasm::parse("qreg q[2];\nx q[5];\n")
        .unwrap_err()
        .to_string();
    assert!(err.contains("out of range"), "{err}");
}

#[test]
fn error_unknown_register() {
    let err = qasm::parse("qreg q[2];\nx r[0];\n")
        .unwrap_err()
        .to_string();
    assert!(err.contains("unknown register `r`"), "{err}");
}

// Regression: malformed bracket order used to panic (slice out of order)
// instead of returning an error.
#[test]
fn error_reversed_brackets_in_qreg() {
    let err = qasm::parse("qreg q]3[;\n").unwrap_err().to_string();
    assert!(err.contains("malformed register declaration"), "{err}");
}

#[test]
fn error_reversed_brackets_in_operand() {
    let err = qasm::parse("qreg q[2];\nx q]0[;\n")
        .unwrap_err()
        .to_string();
    assert!(err.contains("malformed qubit reference"), "{err}");
}

// Regression: a huge register used to abort the process on allocation.
#[test]
fn error_register_too_large() {
    let err = qasm::parse("qreg q[40];\n").unwrap_err().to_string();
    assert!(err.contains("exceeds the maximum"), "{err}");
}

// Regression (found by cargo fuzz): summing multiple huge register sizes used
// to overflow `total` and panic under debug assertions.
#[test]
fn error_register_sizes_do_not_overflow() {
    let huge = usize::MAX;
    let err = qasm::parse(&format!("qreg a[{huge}];\nqreg b[{huge}];\n"))
        .unwrap_err()
        .to_string();
    assert!(err.contains("exceeds the maximum"), "{err}");
}

#[test]
fn bundled_bell_qasm_matches_builder() {
    let src = include_str!("../programs/bell.qasm");
    let imported = qasm::parse(src)
        .expect("bundled example should parse")
        .run();

    let mut expected = Circuit::new(2);
    expected.h(0).cnot(0, 1);
    assert_same_probs(&imported, &expected.run());
}

/// Regression: a `gate` body's `{ ... }` holds its own semicolons and ends
/// without one, so splitting statements naively on `;` glued the closing brace
/// onto the next statement. Since `gate` blocks conventionally come *before*
/// `qreg` — that is how Qiskit emits them — the register declaration was
/// swallowed and the parser reported a missing `qreg` for a file that had one.
#[test]
fn gate_block_before_qreg_finds_the_register() {
    let src = "OPENQASM 2.0;\ngate g a { h a; }\nqreg q[1];\nh q[0];\n";
    let c = qasm::parse(src).expect("the qreg after a gate body must be seen");
    assert_eq!(c.run().n_qubits(), 1);
}

/// Unbalanced braces are reported as such, not as some downstream confusion.
#[test]
fn unbalanced_braces_are_reported() {
    let open = qasm::parse("OPENQASM 2.0;\ngate g a { h a;\nqreg q[1];\n")
        .unwrap_err()
        .to_string();
    assert!(open.contains("unterminated `{`"), "{open}");

    let close = qasm::parse("OPENQASM 2.0;\nqreg q[1];\n}\n")
        .unwrap_err()
        .to_string();
    assert!(close.contains("unexpected `}`"), "{close}");
}

/// A `gate` declaration is expanded at its call site, so a circuit written
/// with one matches the same circuit written with the primitives.
#[test]
fn gate_declaration_expands() {
    let src = "\
OPENQASM 2.0;
gate bell a,b { h a; cx a,b; }
qreg q[2];
bell q[0],q[1];
";
    let expanded = qasm::parse(src).expect("should parse");
    let direct = qasm::parse("qreg q[2];\nh q[0];\ncx q[0],q[1];\n").unwrap();
    assert_same_amplitudes(&expanded.run(), &direct.run());
}

/// Parameters are substituted into the body's angle expressions, which may be
/// arithmetic over them.
#[test]
fn gate_declaration_takes_parameters() {
    let src = "\
OPENQASM 2.0;
gate myrz(theta) a { u1(theta/2) a; u1(theta/2) a; }
qreg q[1];
h q[0];
myrz(0.8) q[0];
";
    let expanded = qasm::parse(src).expect("should parse");
    let direct = qasm::parse("qreg q[1];\nh q[0];\nu1(0.8) q[0];\n").unwrap();
    assert_same_amplitudes(&expanded.run(), &direct.run());
}

/// A body may call another declaration, and qubit arguments are remapped at
/// each level — including when the caller passes them in a different order.
#[test]
fn gate_declarations_nest_and_remap_qubits() {
    let src = "\
OPENQASM 2.0;
gate flip a { x a; }
gate pair a,b { flip b; cx b,a; }
qreg q[2];
pair q[0],q[1];
";
    let expanded = qasm::parse(src).expect("should parse");
    let direct = qasm::parse("qreg q[2];\nx q[1];\ncx q[1],q[0];\n").unwrap();
    assert_same_amplitudes(&expanded.run(), &direct.run());
}

/// The OpenQASM primitives `U` and `CX` work, so a program that defines the
/// standard gates itself — the way `qelib1.inc` does — runs.
#[test]
fn openqasm_primitives_are_supported() {
    let src = "\
OPENQASM 2.0;
gate myh a { U(pi/2,0,pi) a; }
gate mycx c,t { CX c,t; }
qreg q[2];
myh q[0];
mycx q[0],q[1];
";
    let expanded = qasm::parse(src).expect("should parse");
    let direct = qasm::parse("qreg q[2];\nh q[0];\ncx q[0],q[1];\n").unwrap();
    assert_same_amplitudes(&expanded.run(), &direct.run());
}

/// A declaration shadows the built-in of the same name, so the file's own
/// definition is what runs.
#[test]
fn declaration_shadows_the_builtin() {
    // Redefine `x` as a no-op; the qubit must stay in |0>.
    let src = "OPENQASM 2.0;\ngate x a { id a; }\nqreg q[1];\nx q[0];\n";
    let c = qasm::parse(src).expect("should parse");
    assert!(c.run().probability(0) > 1.0 - 1e-12);
}

/// Malformed and abusive declarations are rejected with a message that names
/// the problem — including the two guards against a runaway expansion.
#[test]
fn gate_declaration_errors() {
    let cases = [
        // Wrong number of qubits at the call site.
        (
            "gate g a,b { cx a,b; }\nqreg q[2];\ng q[0];",
            "takes 2 qubit(s) and 0 angle(s)",
        ),
        // Wrong number of angles.
        (
            "gate g(t) a { u1(t) a; }\nqreg q[1];\ng q[0];",
            "takes 1 qubit(s) and 1 angle(s)",
        ),
        // The body names a qubit that is not a formal argument.
        (
            "gate g a { cx a,zz; }\nqreg q[2];\ng q[0];",
            "not a qubit argument",
        ),
        // The body uses an unknown angle name.
        (
            "gate g a { u1(nope) a; }\nqreg q[1];\ng q[0];",
            "unknown name",
        ),
        (
            "gate g a { h a; }\ngate g a { x a; }\nqreg q[1];\ng q[0];",
            "duplicate `gate g`",
        ),
        ("gate g { h a; }\nqreg q[1];", "at least one qubit argument"),
        (
            "gate g a,a { h a; }\nqreg q[1];\ng q[0];",
            "duplicate qubit argument",
        ),
        // Direct recursion must terminate rather than blow the stack.
        ("gate g a { g a; }\nqreg q[1];\ng q[0];", "recursive"),
    ];
    for (src, needle) in cases {
        let err = qasm::parse(src).unwrap_err().to_string();
        assert!(err.contains(needle), "for `{src}`: {err}");
    }
}

/// A short file can describe an exponential expansion, so the budget — not the
/// input length — is what bounds the work.
#[test]
fn exponential_expansion_is_bounded() {
    let mut src = String::from("gate g0 a { x a; }\n");
    for i in 1..40 {
        src.push_str(&format!("gate g{i} a {{ g{} a; g{} a; }}\n", i - 1, i - 1));
    }
    src.push_str("qreg q[1];\ng39 q[0];\n");

    let err = qasm::parse(&src).unwrap_err().to_string();
    assert!(err.contains("expands to more than"), "{err}");
}

/// The gates `qelib1.inc` would supply are implemented natively, since
/// `include` is ignored. Each is checked against Qiskit's unitary separately;
/// here we pin the identities that hold within this engine.
#[test]
fn qelib1_gates_import() {
    // `u` is `u3`, and `u0` is an idle of some duration — an identity here.
    assert_same_amplitudes(
        &qasm::parse("qreg q[1];\nu(0.6,-0.3,1.2) q[0];\n")
            .unwrap()
            .run(),
        &qasm::parse("qreg q[1];\nu3(0.6,-0.3,1.2) q[0];\n")
            .unwrap()
            .run(),
    );
    assert_same_amplitudes(
        &qasm::parse("qreg q[1];\nh q[0];\nu0(3) q[0];\n")
            .unwrap()
            .run(),
        &qasm::parse("qreg q[1];\nh q[0];\n").unwrap().run(),
    );

    // `sx` twice is X, and `sxdg` undoes it.
    let c = qasm::parse("qreg q[1];\nsx q[0];\nsx q[0];\n").unwrap();
    assert_relative_eq!(c.run().probability(1), 1.0, epsilon = 1e-12);
    let c = qasm::parse("qreg q[1];\nh q[0];\nsx q[0];\nsxdg q[0];\n").unwrap();
    assert_same_amplitudes(
        &c.run(),
        &qasm::parse("qreg q[1];\nh q[0];\n").unwrap().run(),
    );

    // `c3x`/`c4x` are true multi-controlled X: they fire only on all-ones.
    let c = qasm::parse(
        "qreg q[5];\nx q[0];\nx q[1];\nx q[2];\nx q[3];\nc4x q[0],q[1],q[2],q[3],q[4];\n",
    )
    .unwrap();
    assert_relative_eq!(c.run().probability(0b11111), 1.0, epsilon = 1e-12);
    let c = qasm::parse("qreg q[5];\nx q[0];\nx q[1];\nc4x q[0],q[1],q[2],q[3],q[4];\n").unwrap();
    assert_relative_eq!(c.run().probability(0b00011), 1.0, epsilon = 1e-12);

    // `c3sqrtx` applied twice is `c3x`.
    let twice = qasm::parse(
        "qreg q[4];\nh q[0];\nh q[1];\nh q[2];\nc3sqrtx q[0],q[1],q[2],q[3];\nc3sqrtx q[0],q[1],q[2],q[3];\n",
    ).unwrap();
    let once =
        qasm::parse("qreg q[4];\nh q[0];\nh q[1];\nh q[2];\nc3x q[0],q[1],q[2],q[3];\n").unwrap();
    assert_same_amplitudes(&twice.run(), &once.run());

    // The relative-phase Toffoli permutes like `ccx` but is not equal to it:
    // it differs by phases, which is the whole point of the cheaper form.
    let rccx = qasm::parse("qreg q[3];\nh q[0];\nh q[1];\nrccx q[0],q[1],q[2];\n").unwrap();
    let ccx = qasm::parse("qreg q[3];\nh q[0];\nh q[1];\nccx q[0],q[1],q[2];\n").unwrap();
    for i in 0..8 {
        assert_relative_eq!(
            rccx.run().probability(i),
            ccx.run().probability(i),
            epsilon = 1e-12
        );
    }
    assert!(
        rccx.run()
            .amplitudes()
            .iter()
            .zip(ccx.run().amplitudes())
            .any(|(a, b)| (a - b).norm() > 1e-9),
        "rccx must differ from ccx in phase"
    );
}

/// `rxx`/`rzz` follow qelib1's decomposition, which differs from Qiskit's gate
/// object by a global phase of theta/2 — unobservable, and expressible neither
/// in OpenQASM 2 nor in this engine's circuit representation.
#[test]
fn two_qubit_rotations_import() {
    // rzz(theta) is diagonal: it leaves populations alone.
    let c = qasm::parse("qreg q[2];\nh q[0];\nh q[1];\nrzz(0.7) q[0],q[1];\n").unwrap();
    for i in 0..4 {
        assert_relative_eq!(c.run().probability(i), 0.25, epsilon = 1e-12);
    }
    // rzz(0) and rxx(0) are identities.
    for src in ["rzz(0) q[0],q[1];", "rxx(0) q[0],q[1];"] {
        let c = qasm::parse(&format!("qreg q[2];\nh q[0];\n{src}\n")).unwrap();
        assert_same_amplitudes(
            &c.run(),
            &qasm::parse("qreg q[2];\nh q[0];\n").unwrap().run(),
        );
    }
}

/// `measure` is honoured rather than ignored: it collapses the qubit, so a
/// gate after it sees a definite state. Previously this whole circuit reduced
/// to `h; h` and reported |0> with certainty.
#[test]
fn measure_collapses_on_import() {
    let c = qasm::parse(
        "OPENQASM 2.0;\nqreg q[1];\ncreg c[1];\nh q[0];\nmeasure q[0] -> c[0];\nh q[0];\n",
    )
    .expect("should parse");
    let hist = c.sample(2000, 5);
    let zeros = hist.get(&0).copied().unwrap_or(0);
    assert!(
        (800..=1200).contains(&zeros),
        "expected a coin flip, got {zeros}/2000 zeros"
    );
}

/// Malformed measurements are rejected rather than half-understood.
#[test]
fn measure_errors() {
    for (src, needle) in [
        (
            "qreg q[1];\ncreg c[1];\nmeasure q[0];\n",
            "needs `-> target`",
        ),
        (
            "qreg q[1];\ncreg c[1];\nmeasure q[0] ->;\n",
            "classical target",
        ),
        (
            "qreg q[2];\ncreg c[2];\nmeasure q[0],q[1] -> c[0];\n",
            "expected one qubit",
        ),
    ] {
        let err = qasm::parse(src).unwrap_err().to_string();
        assert!(err.contains(needle), "for `{src}`: {err}");
    }
}

/// `reset` is honoured on import, in both the indexed and whole-register forms.
#[test]
fn reset_imports() {
    let indexed = qasm::parse("qreg q[2];\nh q[0];\nh q[1];\nreset q[0];\n").expect("parse");
    let state = indexed.run();
    assert_relative_eq!(state.prob_qubit_one(0), 0.0, epsilon = 1e-12);

    let whole = qasm::parse("qreg q[2];\nh q[0];\nh q[1];\nreset q;\n").expect("parse");
    assert_relative_eq!(whole.run().probability(0b00), 1.0, epsilon = 1e-12);

    let err = qasm::parse("qreg q[1];\nreset nope;\n")
        .unwrap_err()
        .to_string();
    assert!(err.contains("unknown register"), "{err}");
}

/// `if (c == value)` guards a gate on the classical register, so a measurement
/// can steer what happens next.
#[test]
fn conditional_imports() {
    let src = "\
OPENQASM 2.0;
qreg q[2];
creg c[2];
h q[0];
measure q[0] -> c[0];
if(c==1) x q[1];
measure q[1] -> c[1];
";
    let c = qasm::parse(src).expect("should parse");
    let hist = c.sample(2000, 4);
    // The conditional flip correlates the qubits: only 00 and 11 occur.
    assert_eq!(hist.get(&0b01).copied().unwrap_or(0), 0);
    assert_eq!(hist.get(&0b10).copied().unwrap_or(0), 0);
}

/// A measurement may name a classical bit other than its own index.
#[test]
fn measure_honours_the_classical_target() {
    let src = "qreg q[2];\ncreg c[2];\nx q[1];\nmeasure q[1] -> c[0];\nif(c==1) x q[0];\n";
    let state = qasm::parse(src).expect("should parse").run();
    assert_relative_eq!(state.probability(0b11), 1.0, epsilon = 1e-12);
}

/// Malformed conditionals and unsupported shapes are named, not guessed at.
#[test]
fn conditional_errors() {
    for (src, needle) in [
        (
            "qreg q[1];\ncreg c[1];\nif c==0 x q[0];\n",
            "needs `(c==value)`",
        ),
        (
            "qreg q[1];\ncreg c[1];\nif(c=0) x q[0];\n",
            "compare with `==`",
        ),
        (
            "qreg q[1];\ncreg c[1];\nif(d==0) x q[0];\n",
            "unknown classical register",
        ),
        (
            "qreg q[1];\ncreg c[1];\nif(c==x) x q[0];\n",
            "non-negative integer",
        ),
        ("qreg q[1];\ncreg c[1];\nif(c==0)\n", "no statement"),
        (
            "qreg q[1];\ncreg c[1];\nif(c==0) measure q[0] -> c[0];\n",
            "only a gate may be conditional",
        ),
        (
            "qreg q[1];\ncreg c[4];\n",
            "more than the 1 qubits this engine models",
        ),
        (
            "qreg q[2];\ncreg c[2];\nmeasure q -> c[0];\n",
            "onto 1 classical bit",
        ),
    ] {
        let err = qasm::parse(src).unwrap_err().to_string();
        assert!(err.contains(needle), "for `{src}`: {err}");
    }
}

/// The controlled qelib1 gates, checked against the engine's own primitives.
/// These were reached only by the cross-validation harness, which runs the
/// compiled binary in a subprocess and so contributes no in-tree coverage: a
/// refactor could break one and `cargo test` would stay green.
#[test]
fn controlled_qelib1_gates_match_their_primitives() {
    // crx/cry with the control set are the plain rotation on the target.
    for (name, builder) in [("crx", "rx"), ("cry", "ry")] {
        let controlled = qasm::parse(&format!("qreg q[2];\nx q[0];\n{name}(0.7) q[0],q[1];\n"))
            .expect("should parse");
        let mut plain = Circuit::new(2);
        plain.x(0);
        match builder {
            "rx" => plain.rx(0.7, 1),
            _ => plain.ry(0.7, 1),
        };
        assert_same_amplitudes(&controlled.run(), &plain.run());

        // With the control clear they do nothing.
        let idle = qasm::parse(&format!("qreg q[2];\n{name}(0.7) q[0],q[1];\n")).unwrap();
        assert_relative_eq!(idle.run().probability(0), 1.0, epsilon = 1e-12);
    }

    // csx twice with the control set is a CNOT, since sx squared is X.
    let twice = qasm::parse("qreg q[2];\nx q[0];\ncsx q[0],q[1];\ncsx q[0],q[1];\n").unwrap();
    let mut cnot = Circuit::new(2);
    cnot.x(0).cnot(0, 1);
    assert_same_amplitudes(&twice.run(), &cnot.run());

    // cu(theta,phi,lambda,gamma) is a controlled-u3 plus a phase on the control.
    let cu = qasm::parse("qreg q[2];\nh q[0];\ncu(0.6,-0.3,1.2,0.4) q[0],q[1];\n").unwrap();
    let mut equivalent = Circuit::new(2);
    equivalent.h(0).p(0.4, 0).cu3(0.6, -0.3, 1.2, 0, 1);
    assert_same_amplitudes(&cu.run(), &equivalent.run());
}

/// `rc3x` is a *relative-phase* three-control X: it permutes like `c3x` but
/// differs in phase, which is the point of the cheaper form.
#[test]
fn rc3x_matches_c3x_in_populations_but_not_phase() {
    let prep = "qreg q[4];\nh q[0];\nh q[1];\nh q[2];\n";
    let relative = qasm::parse(&format!("{prep}rc3x q[0],q[1],q[2],q[3];\n")).unwrap();
    let exact = qasm::parse(&format!("{prep}c3x q[0],q[1],q[2],q[3];\n")).unwrap();

    let (a, b) = (relative.run(), exact.run());
    for i in 0..16 {
        assert_relative_eq!(a.probability(i), b.probability(i), epsilon = 1e-12);
    }
    assert!(
        a.amplitudes()
            .iter()
            .zip(b.amplitudes())
            .any(|(x, y)| (x - y).norm() > 1e-9),
        "rc3x must differ from c3x in phase"
    );
}

/// Declaration errors that no valid program reaches.
#[test]
fn declaration_errors() {
    for (src, needle) in [
        ("qreg q[1];\nqreg q[1];\n", "duplicate register"),
        (
            "qreg q[1];\ncreg c[1];\ncreg c[1];\n",
            "duplicate classical register",
        ),
        ("qreg q[0];\n", "size must be >= 1"),
        ("qreg q[1];\ngate { h a; }\n", "needs a name"),
        ("qreg q[1];\ngate g(1) a { h a; }\n", "invalid parameter"),
        ("qreg q[2];\ncx q[0],q[0];\n", "must be distinct qubits"),
        ("qreg q[1];\nreset ;\n", "expected a qubit"),
    ] {
        let err = qasm::parse(src).unwrap_err().to_string();
        assert!(err.contains(needle), "for `{src}`: {err}");
    }
}
