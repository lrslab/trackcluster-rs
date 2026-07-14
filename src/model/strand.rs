use thiserror::Error;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
/// Genomic strand parsed from BED's `+`, `-`, or `.` values.
pub enum Strand {
    /// Forward (`+`) strand.
    Plus,
    /// Reverse (`-`) strand.
    Minus,
    /// Unknown or unstranded (`.`).
    Unknown,
}

#[derive(Error, Debug)]
/// Error returned for a non-BED strand token.
pub enum StrandParseError {
    /// The supplied token was not `+`, `-`, or `.`.
    #[error("invalid strand {value:?}")]
    Invalid {
        /// Supplied strand token.
        value: String,
    },
}

impl Strand {
    /// Return the canonical BED character for this strand.
    pub fn as_char(self) -> char {
        match self {
            Self::Plus => '+',
            Self::Minus => '-',
            Self::Unknown => '.',
        }
    }
}

impl TryFrom<&str> for Strand {
    type Error = StrandParseError;

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        match value {
            "+" => Ok(Self::Plus),
            "-" => Ok(Self::Minus),
            "." => Ok(Self::Unknown),
            _ => Err(StrandParseError::Invalid {
                value: value.to_owned(),
            }),
        }
    }
}
