//! Validated domain types for scientific and runtime configuration.

use std::num::{NonZeroU32, NonZeroUsize};

use thiserror::Error;

/// Error returned when a scientific or runtime parameter is outside its domain.
#[derive(Clone, Debug, Error, PartialEq)]
pub enum ParameterError {
    /// An input identity violates the clustering catalog contract.
    #[error("invalid input identity: {message}")]
    InvalidIdentity {
        /// Detailed identity validation failure.
        message: String,
    },
    /// A floating-point value is NaN or infinite.
    #[error("{parameter} must be finite; got {value}")]
    NonFinite {
        /// Parameter name shown to the user.
        parameter: &'static str,
        /// Rejected value.
        value: f64,
    },
    /// A fraction or distance cutoff is outside the closed unit interval.
    #[error("{parameter} must be within [0, 1]; got {value}")]
    OutsideUnitInterval {
        /// Parameter name shown to the user.
        parameter: &'static str,
        /// Rejected value.
        value: f64,
    },
    /// A finite numeric value is negative.
    #[error("{parameter} must be nonnegative; got {value}")]
    Negative {
        /// Parameter name shown to the user.
        parameter: &'static str,
        /// Rejected value.
        value: f64,
    },
    /// A count that must be positive was zero.
    #[error("{parameter} must be at least 1; got 0")]
    Zero {
        /// Parameter name shown to the user.
        parameter: &'static str,
    },
}

impl ParameterError {
    pub(crate) fn invalid_identity(error: impl std::fmt::Display) -> Self {
        Self::InvalidIdentity {
            message: error.to_string(),
        }
    }
}

/// A finite fraction or normalized distance cutoff in `[0, 1]`.
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct UnitFraction(f64);

impl UnitFraction {
    /// Validate and construct a unit fraction.
    pub fn new(parameter: &'static str, value: f64) -> Result<Self, ParameterError> {
        if !value.is_finite() {
            return Err(ParameterError::NonFinite { parameter, value });
        }
        if !(0.0..=1.0).contains(&value) {
            return Err(ParameterError::OutsideUnitInterval { parameter, value });
        }
        Ok(Self(value))
    }

    /// Return the validated primitive value.
    pub const fn get(self) -> f64 {
        self.0
    }
}

/// A finite, nonnegative floating-point weight.
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct NonNegativeWeight(f64);

impl NonNegativeWeight {
    /// Validate and construct a nonnegative weight.
    pub fn new(parameter: &'static str, value: f64) -> Result<Self, ParameterError> {
        if !value.is_finite() {
            return Err(ParameterError::NonFinite { parameter, value });
        }
        if value < 0.0 {
            return Err(ParameterError::Negative { parameter, value });
        }
        Ok(Self(value))
    }

    /// Return the validated primitive value.
    pub const fn get(self) -> f64 {
        self.0
    }
}

/// A validated, nonzero worker count.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct WorkerThreads(NonZeroUsize);

impl WorkerThreads {
    /// Validate and construct a worker count.
    pub fn new(value: usize) -> Result<Self, ParameterError> {
        NonZeroUsize::new(value)
            .map(Self)
            .ok_or(ParameterError::Zero {
                parameter: "worker threads",
            })
    }

    /// Return the validated worker count.
    pub const fn get(self) -> usize {
        self.0.get()
    }
}

/// A validated, nonzero number of iterative clustering rounds.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct BatchRounds(NonZeroUsize);

impl BatchRounds {
    /// Validate and construct a batch-round limit.
    pub fn new(value: usize) -> Result<Self, ParameterError> {
        NonZeroUsize::new(value)
            .map(Self)
            .ok_or(ParameterError::Zero {
                parameter: "batch rounds",
            })
    }

    /// Return the validated round limit.
    pub const fn get(self) -> usize {
        self.0.get()
    }
}

/// A validated base-pair offset.
///
/// The contained `u32` makes negative and non-integral offsets unrepresentable.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct BasePairOffset(u32);

impl BasePairOffset {
    /// Construct a nonnegative base-pair offset.
    pub const fn new(value: u32) -> Self {
        Self(value)
    }

    /// Return the offset in base pairs.
    pub const fn get(self) -> u32 {
        self.0
    }
}

/// A nonzero read-support threshold.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct MinimumSupport(NonZeroUsize);

impl MinimumSupport {
    /// Validate and construct a read-support threshold.
    pub fn new(parameter: &'static str, value: usize) -> Result<Self, ParameterError> {
        NonZeroUsize::new(value)
            .map(Self)
            .ok_or(ParameterError::Zero { parameter })
    }

    /// Return the support threshold.
    pub const fn get(self) -> usize {
        self.0.get()
    }
}

/// A nonzero weighted-support threshold.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct WeightedMinimumSupport(NonZeroU32);

impl WeightedMinimumSupport {
    /// Validate and construct a weighted-support threshold.
    pub fn new(parameter: &'static str, value: u32) -> Result<Self, ParameterError> {
        NonZeroU32::new(value)
            .map(Self)
            .ok_or(ParameterError::Zero { parameter })
    }

    /// Return the support threshold.
    pub const fn get(self) -> u32 {
        self.0.get()
    }
}

fn parse_f64(value: &str, parameter: &'static str) -> Result<f64, String> {
    value
        .parse::<f64>()
        .map_err(|error| format!("{parameter} must be a number: {error}"))
}

/// Clap value parser for fractions and normalized distance cutoffs.
pub fn parse_unit_fraction(value: &str) -> Result<f64, String> {
    let parsed = parse_f64(value, "fraction or distance cutoff")?;
    UnitFraction::new("fraction or distance cutoff", parsed)
        .map(UnitFraction::get)
        .map_err(|error| error.to_string())
}

/// Clap value parser for the overlap intron weight.
pub fn parse_nonnegative_weight(value: &str) -> Result<f64, String> {
    let parsed = parse_f64(value, "intron weight")?;
    NonNegativeWeight::new("intron weight", parsed)
        .map(NonNegativeWeight::get)
        .map_err(|error| error.to_string())
}

/// Clap value parser for worker threads.
pub fn parse_worker_threads(value: &str) -> Result<usize, String> {
    let parsed = value
        .parse::<usize>()
        .map_err(|error| format!("worker threads must be a positive integer: {error}"))?;
    WorkerThreads::new(parsed)
        .map(WorkerThreads::get)
        .map_err(|error| error.to_string())
}

/// Clap value parser for iterative clustering rounds.
pub fn parse_batch_rounds(value: &str) -> Result<usize, String> {
    let parsed = value
        .parse::<usize>()
        .map_err(|error| format!("batch rounds must be a positive integer: {error}"))?;
    BatchRounds::new(parsed)
        .map(BatchRounds::get)
        .map_err(|error| error.to_string())
}

/// Clap value parser for a read-support threshold.
pub fn parse_minimum_support(value: &str) -> Result<usize, String> {
    let parsed = value
        .parse::<usize>()
        .map_err(|error| format!("minimum support must be a positive integer: {error}"))?;
    MinimumSupport::new("minimum support", parsed)
        .map(MinimumSupport::get)
        .map_err(|error| error.to_string())
}

/// Clap value parser for a weighted-support threshold.
pub fn parse_weighted_minimum_support(value: &str) -> Result<u32, String> {
    let parsed = value
        .parse::<u32>()
        .map_err(|error| format!("minimum support must be a positive integer: {error}"))?;
    WeightedMinimumSupport::new("minimum support", parsed)
        .map(WeightedMinimumSupport::get)
        .map_err(|error| error.to_string())
}

/// Clap value parser for base-pair offsets.
pub fn parse_base_pair_offset(value: &str) -> Result<u32, String> {
    value
        .parse::<u32>()
        .map(BasePairOffset::new)
        .map(BasePairOffset::get)
        .map_err(|error| format!("base-pair offset must be a nonnegative integer: {error}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unit_fraction_accepts_closed_interval_only() {
        assert_eq!(UnitFraction::new("cutoff", 0.0).unwrap().get(), 0.0);
        assert_eq!(UnitFraction::new("cutoff", 1.0).unwrap().get(), 1.0);
        for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -0.01, 1.01] {
            assert!(UnitFraction::new("cutoff", invalid).is_err(), "{invalid}");
        }
    }

    #[test]
    fn nonnegative_weight_requires_finite_value() {
        assert_eq!(NonNegativeWeight::new("weight", 0.0).unwrap().get(), 0.0);
        for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -0.01] {
            assert!(
                NonNegativeWeight::new("weight", invalid).is_err(),
                "{invalid}"
            );
        }
    }

    #[test]
    fn nonzero_count_types_reject_zero() {
        assert!(WorkerThreads::new(0).is_err());
        assert!(BatchRounds::new(0).is_err());
        assert!(MinimumSupport::new("support", 0).is_err());
        assert!(WeightedMinimumSupport::new("support", 0).is_err());
        assert_eq!(WorkerThreads::new(4).unwrap().get(), 4);
        assert_eq!(BatchRounds::new(7).unwrap().get(), 7);
        assert_eq!(MinimumSupport::new("support", 2).unwrap().get(), 2);
        assert_eq!(WeightedMinimumSupport::new("support", 3).unwrap().get(), 3);
    }

    #[test]
    fn clap_parsers_reject_nonfinite_and_zero_values() {
        for invalid in ["NaN", "inf", "-inf", "-0.1", "1.1"] {
            assert!(parse_unit_fraction(invalid).is_err(), "{invalid}");
        }
        for invalid in ["NaN", "inf", "-inf", "-0.1"] {
            assert!(parse_nonnegative_weight(invalid).is_err(), "{invalid}");
        }
        assert!(parse_worker_threads("0").is_err());
        assert!(parse_batch_rounds("0").is_err());
        assert!(parse_minimum_support("0").is_err());
        assert!(parse_weighted_minimum_support("0").is_err());
    }
}
