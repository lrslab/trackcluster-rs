//! Strict genomics file readers and writers.

/// Pure-Rust BAM to TrackCluster bigGenePred conversion.
pub mod bam;
/// Deterministic BAM splitting by query-name membership.
pub mod bam_subset;
/// BED12 and bigGenePred-compatible parsing and serialization.
pub mod bed;
/// Exact genomic-base coverage from primary BAM alignments.
pub mod coverage;
/// Indexed FASTA reference-base access.
pub mod fasta;
/// GFF3/GTF transcript annotation parsing and bigGenePred conversion.
pub mod gff;
/// GTF, GFF3, and SQANTI3-oriented transcript exports.
pub mod interchange;
/// m6Anet RNA002 read-probability import and transcript-to-genome projection.
pub mod m6anet;
/// Tab-delimited multi-sample manifest parsing.
pub mod manifest;
/// Normalized modification observation and assay metadata I/O.
pub mod mod_calls;
/// Modification sample manifest parsing.
pub mod mod_manifest;
/// Dorado/modBAM MM/ML decoding and query-to-genome projection.
pub mod modbam;
