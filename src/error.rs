//! Error types for the text and OpenQASM front ends.

use std::fmt;

/// An error from parsing a text program ([`crate::program::parse`]) or
/// OpenQASM source ([`crate::qasm::parse`]).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParseError {
    /// 1-based source line, when the parser tracks one (the text program
    /// format does; the OpenQASM importer reports the offending statement in
    /// the message instead).
    pub line: Option<usize>,
    /// Human-readable description of the problem.
    pub message: String,
}

impl ParseError {
    /// A parse error with no associated line.
    pub(crate) fn new(message: impl Into<String>) -> Self {
        ParseError {
            line: None,
            message: message.into(),
        }
    }

    /// A parse error located at a 1-based source line.
    pub(crate) fn at_line(line: usize, message: impl Into<String>) -> Self {
        ParseError {
            line: Some(line),
            message: message.into(),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.line {
            Some(line) => write!(f, "line {line}: {}", self.message),
            None => f.write_str(&self.message),
        }
    }
}

impl std::error::Error for ParseError {}
