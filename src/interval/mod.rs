pub mod cluster_span;
pub mod index;
#[cfg(feature = "index-binned")]
pub mod index_binned;
pub mod intersect_sweep;
pub mod partition;
pub mod refine;
pub mod sort;

use crate::model::Strand;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum StrandMode {
    #[default]
    Ignore,
    Match,
}

impl StrandMode {
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
