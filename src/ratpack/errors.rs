// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

//! Rust translation of `CalcErr.h`.
//! Error codes thrown by ratpak and caught by Calculator.

use std::fmt;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CalcError {
    /// The current operation would require a divide by zero to complete.
    DivideByZero,
    /// The given input is not within the domain of this function.
    Domain,
    /// The result of this function is undefined.
    Indefinite,
    /// The result of this function is Positive Infinity.
    PosInfinity,
    /// The result of this function is Negative Infinity.
    NegInfinity,
    /// The given input is within the domain but beyond calc's computable range.
    InvalidRange,
    /// There is not enough free memory to complete the requested function.
    OutOfMemory,
    /// The result of this operation is an overflow.
    Overflow,
    /// The result of this operation is undefined.
    NoResult,
}

impl fmt::Display for CalcError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CalcError::DivideByZero => write!(f, "Divide by zero"),
            CalcError::Domain => write!(f, "Input out of domain"),
            CalcError::Indefinite => write!(f, "Indefinite result"),
            CalcError::PosInfinity => write!(f, "Positive infinity"),
            CalcError::NegInfinity => write!(f, "Negative infinity"),
            CalcError::InvalidRange => write!(f, "Invalid range"),
            CalcError::OutOfMemory => write!(f, "Out of memory"),
            CalcError::Overflow => write!(f, "Overflow"),
            CalcError::NoResult => write!(f, "No result"),
        }
    }
}

impl std::error::Error for CalcError {}

pub type CalcResult<T> = Result<T, CalcError>;
