#![no_main]

//! Fuzz the text program parser: for any input, `program::parse` must return a
//! `Result` and never panic.

use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    if let Ok(text) = std::str::from_utf8(data) {
        let _ = qsimulator::program::parse(text);
    }
});
