//! Strict genomics file readers and writers.

/// Pure-Rust BAM to TrackCluster bigGenePred conversion.
pub mod bam;
/// BED12 and bigGenePred-compatible parsing and serialization.
pub mod bed;
/// GFF3/GTF transcript annotation parsing and bigGenePred conversion.
pub mod gff;
/// GTF, GFF3, and SQANTI3-oriented transcript exports.
pub mod interchange;
/// Tab-delimited multi-sample manifest parsing.
pub mod manifest;
