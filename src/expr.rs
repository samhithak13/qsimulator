//! Angle expression evaluation, shared by both front ends.
//!
//! OpenQASM 2 angles are arithmetic expressions, not just literals: a `gate`
//! body writes things like `rz(-theta/2)` and `p(pi/8)` over the declaration's
//! formal parameters. This evaluates that grammar against a variable
//! environment.
//!
//! Supported: decimal and exponent-form numbers, `pi`, variables, the unary and
//! binary operators `+ - * / ^` with the usual precedence (`^` binds tightest
//! and associates to the right), parentheses, and the OpenQASM functions `sin`,
//! `cos`, `tan`, `exp`, `ln`, `sqrt`.
//!
//! For compatibility with the native text format, a number directly against a
//! name or parenthesis is an implicit multiplication, so `2pi` means `2*pi`.

use std::collections::HashMap;
use std::f64::consts::PI;

/// Evaluate an angle expression, resolving names through `vars`.
///
/// `vars` is empty for a top-level angle and holds the formal parameters when
/// expanding the body of a `gate` declaration.
pub(crate) fn eval(src: &str, vars: &HashMap<String, f64>) -> Result<f64, String> {
    let mut p = Parser {
        chars: src.char_indices().peekable(),
        src,
        vars,
    };
    let value = p.expr()?;
    p.skip_space();
    if let Some((_, c)) = p.chars.peek() {
        return Err(format!("invalid angle `{src}` (unexpected `{c}`)"));
    }
    if !value.is_finite() {
        return Err(format!("invalid angle `{src}` (not a finite number)"));
    }
    Ok(value)
}

struct Parser<'a> {
    chars: std::iter::Peekable<std::str::CharIndices<'a>>,
    src: &'a str,
    vars: &'a HashMap<String, f64>,
}

impl Parser<'_> {
    fn invalid(&self) -> String {
        format!("invalid angle `{}`", self.src)
    }

    fn skip_space(&mut self) {
        while matches!(self.chars.peek(), Some((_, c)) if c.is_whitespace()) {
            self.chars.next();
        }
    }

    /// Consume `c` if it is next, reporting whether it was there.
    fn eat(&mut self, c: char) -> bool {
        self.skip_space();
        if matches!(self.chars.peek(), Some((_, got)) if *got == c) {
            self.chars.next();
            true
        } else {
            false
        }
    }

    fn peek(&mut self) -> Option<char> {
        self.skip_space();
        self.chars.peek().map(|(_, c)| *c)
    }

    /// `term (('+' | '-') term)*`
    fn expr(&mut self) -> Result<f64, String> {
        let mut value = self.term()?;
        loop {
            if self.eat('+') {
                value += self.term()?;
            } else if self.eat('-') {
                value -= self.term()?;
            } else {
                return Ok(value);
            }
        }
    }

    /// `unary (('*' | '/') unary)*`
    fn term(&mut self) -> Result<f64, String> {
        let mut value = self.unary()?;
        loop {
            if self.eat('*') {
                value *= self.unary()?;
            } else if self.eat('/') {
                let divisor = self.unary()?;
                if divisor == 0.0 {
                    return Err(format!("invalid angle `{}` (division by zero)", self.src));
                }
                value /= divisor;
            } else {
                return Ok(value);
            }
        }
    }

    /// `('-' | '+') unary | power`
    ///
    /// Sits *above* `^`, so `-2^2` is `-(2^2)` as in ordinary mathematics.
    fn unary(&mut self) -> Result<f64, String> {
        if self.eat('-') {
            return Ok(-self.unary()?);
        }
        if self.eat('+') {
            return self.unary();
        }
        self.power()
    }

    /// `atom ('^' unary)?` — right-associative, so `2^3^2` is `2^9`, and the
    /// exponent may carry its own sign (`2^-1`).
    fn power(&mut self) -> Result<f64, String> {
        let base = self.atom()?;
        if self.eat('^') {
            return Ok(base.powf(self.unary()?));
        }
        Ok(base)
    }

    fn atom(&mut self) -> Result<f64, String> {
        match self.peek().ok_or_else(|| self.invalid())? {
            '(' => {
                self.chars.next();
                let value = self.expr()?;
                if !self.eat(')') {
                    return Err(format!("invalid angle `{}` (missing `)`)", self.src));
                }
                Ok(value)
            }
            c if c.is_ascii_digit() || c == '.' => {
                let number = self.number()?;
                // `2pi` and `2(...)` are implicit multiplications, the form the
                // native text program format has always accepted.
                match self.peek() {
                    Some(next) if next.is_ascii_alphabetic() || next == '(' => {
                        Ok(number * self.atom()?)
                    }
                    _ => Ok(number),
                }
            }
            c if c.is_ascii_alphabetic() || c == '_' => self.name(),
            _ => Err(self.invalid()),
        }
    }

    fn number(&mut self) -> Result<f64, String> {
        let mut text = String::new();
        while let Some((_, c)) = self.chars.peek() {
            let c = *c;
            let is_exponent_sign =
                (c == '+' || c == '-') && matches!(text.chars().last(), Some('e') | Some('E'));
            if c.is_ascii_digit() || c == '.' || c == 'e' || c == 'E' || is_exponent_sign {
                // A trailing `e` may start an identifier (`2exp(1)`) rather than
                // an exponent, so only take it when a digit or sign follows.
                if (c == 'e' || c == 'E') && !self.exponent_follows() {
                    break;
                }
                text.push(c);
                self.chars.next();
            } else {
                break;
            }
        }
        text.parse().map_err(|_| self.invalid())
    }

    /// Look past a pending `e`/`E` for the digits or sign that would make it an
    /// exponent.
    fn exponent_follows(&mut self) -> bool {
        let mut lookahead = self.chars.clone();
        lookahead.next(); // the `e`
        match lookahead.next() {
            Some((_, c)) if c.is_ascii_digit() => true,
            Some((_, '+')) | Some((_, '-')) => {
                matches!(lookahead.next(), Some((_, c)) if c.is_ascii_digit())
            }
            _ => false,
        }
    }

    /// A constant, a variable, or a function call.
    fn name(&mut self) -> Result<f64, String> {
        let mut text = String::new();
        while let Some((_, c)) = self.chars.peek() {
            if c.is_ascii_alphanumeric() || *c == '_' {
                text.push(*c);
                self.chars.next();
            } else {
                break;
            }
        }
        if text == "pi" {
            return Ok(PI);
        }
        if let Some(value) = self.vars.get(&text) {
            return Ok(*value);
        }

        // Anything else must be a function call, or it is an unknown name.
        let func = match text.as_str() {
            "sin" | "cos" | "tan" | "exp" | "ln" | "sqrt" => text,
            _ => {
                return Err(format!(
                    "invalid angle `{}` (unknown name `{text}`)",
                    self.src
                ))
            }
        };
        if !self.eat('(') {
            return Err(format!("invalid angle `{}` (`{func}` needs `(`)", self.src));
        }
        let arg = self.expr()?;
        if !self.eat(')') {
            return Err(format!("invalid angle `{}` (missing `)`)", self.src));
        }
        Ok(match func.as_str() {
            "sin" => arg.sin(),
            "cos" => arg.cos(),
            "tan" => arg.tan(),
            "exp" => arg.exp(),
            "ln" => arg.ln(),
            _ => arg.sqrt(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ev(s: &str) -> Result<f64, String> {
        eval(s, &HashMap::new())
    }

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-12
    }

    /// Every form the native text format accepted before expressions existed
    /// must still evaluate the same way.
    #[test]
    fn legacy_angle_forms_still_parse() {
        for (src, want) in [
            ("0.7", 0.7),
            ("-1.5", -1.5),
            ("1e-5", 1e-5),
            ("2E3", 2000.0),
            ("pi", PI),
            ("pi/2", PI / 2.0),
            ("-pi/4", -PI / 4.0),
            ("2pi", 2.0 * PI),
            ("2*pi", 2.0 * PI),
            ("0.5*pi", 0.5 * PI),
        ] {
            assert!(close(ev(src).unwrap(), want), "{src}");
        }
    }

    #[test]
    fn arithmetic_and_precedence() {
        for (src, want) in [
            ("1+2*3", 7.0),
            ("(1+2)*3", 9.0),
            ("-2^2", -4.0),   // `^` binds tighter than unary minus
            ("2^-1", 0.5),    // the exponent may carry its own sign
            ("2^3^2", 512.0), // right-associative
            ("pi/8", PI / 8.0),
            ("pi/2 + pi/2", PI),
            ("--3", 3.0),
            ("+4", 4.0),
            ("  1  +  1  ", 2.0),
        ] {
            assert!(close(ev(src).unwrap(), want), "{src} -> {:?}", ev(src));
        }
    }

    #[test]
    fn functions_evaluate() {
        assert!(close(ev("sin(0)").unwrap(), 0.0));
        assert!(close(ev("cos(0)").unwrap(), 1.0));
        assert!(close(ev("tan(0)").unwrap(), 0.0));
        assert!(close(ev("ln(1)").unwrap(), 0.0));
        assert!(close(ev("exp(0)").unwrap(), 1.0));
        assert!(close(ev("sqrt(4)").unwrap(), 2.0));
        assert!(close(ev("2exp(0)").unwrap(), 2.0)); // `e` starts a name here
    }

    #[test]
    fn variables_resolve() {
        let vars = HashMap::from([("theta".to_string(), 0.8), ("phi".to_string(), -0.2)]);
        assert!(close(eval("theta", &vars).unwrap(), 0.8));
        assert!(close(eval("-theta/2", &vars).unwrap(), -0.4));
        assert!(close(eval("theta+phi", &vars).unwrap(), 0.6));
        assert!(close(eval("2*theta", &vars).unwrap(), 1.6));
    }

    #[test]
    fn errors_name_the_problem() {
        for (src, needle) in [
            ("bogus", "unknown name"),
            ("pi/0", "division by zero"),
            ("(1+2", "missing `)`"),
            ("1+", "invalid angle"),
            ("", "invalid angle"),
            ("1 2", "unexpected"),
            ("sin 1", "needs `(`"),
            ("1/0", "division by zero"),
        ] {
            let err = ev(src).unwrap_err();
            assert!(err.contains(needle), "for `{src}`: {err}");
        }
    }

    /// A variable shadows nothing and an unknown name is never silently zero.
    #[test]
    fn unknown_variable_is_an_error_not_zero() {
        let vars = HashMap::from([("theta".to_string(), 0.8)]);
        assert!(eval("lambda", &vars).is_err());
    }
}
