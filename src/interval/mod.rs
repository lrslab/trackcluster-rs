//! Interval algorithms whose indices always refer to the caller's original slices.

/// Span-based locus clustering.
pub mod cluster_span;
/// Selectable intersection backends.
pub mod index;
#[cfg(feature = "index-binned")]
/// Reusable fixed-bin interval index enabled by the `index-binned` feature.
pub mod index_binned;
/// Sweep-line span intersection.
pub mod intersect_sweep;
/// Chromosome/strand partitioning.
pub mod partition;
/// Exon and junction comparisons.
pub mod refine;
/// Stable genomic sorting.
pub mod sort;

use crate::model::Strand;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
/// Controls whether interval operations partition records by strand.
pub enum StrandMode {
    /// Ignore strand during matching and partitioning.
    #[default]
    Ignore,
    /// Require equal strand values, including equality between unknown strands.
    Match,
}

impl StrandMode {
    /// Return the strand component used in a partition key.
    pub fn key_strand(self, strand: Strand) -> Option<Strand> {
        match self {
            Self::Ignore => None,
            Self::Match => Some(strand),
        }
    }
}

pub use cluster_span::{cluster_by_span, RangeCluster};
pub use index::intersect_pairs;
#[cfg(feature = "index-binned")]
pub use index::intersect_pairs_binned;
pub use intersect_sweep::{sweep_intersect_pairs, IntersectOpts};
pub use partition::{partition, PartitionKey};
pub use refine::{exonic_overlap_bp, junctions_equal, junctions_subset};
pub use sort::sort_by_coord;
