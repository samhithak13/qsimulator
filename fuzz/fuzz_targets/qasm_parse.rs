#![no_main]

//! Fuzz the OpenQASM importer: for any input, `qasm::parse` must return a
//! `Result` and never panic. Only parsing is exercised, so no state vector is
//! allocated regardless of the declared qubit count.

use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    if let Ok(text) = std::str::from_utf8(data) {
        let _ = qsimulator::qasm::parse(text);
    }
});
