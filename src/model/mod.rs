//! Core validated data structures used throughout TrackCluster-rs.

/// Zero-based genomic coordinate type.
pub mod coord;
/// Half-open genomic interval type.
pub mod interval;
/// Typed TrackCluster/bigGenePred metadata accessors and codec.
pub mod metadata;
/// Genomic strand type.
pub mod strand;
/// BED12-compatible transcript model.
pub mod transcript;

pub use coord::Coord;
pub use interval::Interval;
pub use metadata::{BigGenePredAttrs, TrackMetadataMut, TrackMetadataRef, TranscriptGeometry};
pub use strand::Strand;
pub use transcript::{Bed12Attrs, JunctionSignature, Transcript, TranscriptError};
