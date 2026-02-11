use crate::model::Transcript;

use super::IntersectOpts;

/// An intersect backend that produces `(a_index, b_index)` pairs.
pub trait IntervalIndex {
    fn intersect_pairs(
        a: &[Transcript],
        b: &[Transcript],
        opts: &IntersectOpts,
    ) -> Vec<(usize, usize)>;
}

pub struct SweepIndex;

impl IntervalIndex for SweepIndex {
    fn intersect_pairs(
        a: &[Transcript],
        b: &[Transcript],
        opts: &IntersectOpts,
    ) -> Vec<(usize, usize)> {
        super::intersect_sweep::sweep_intersect_pairs(a, b, opts)
    }
}

#[cfg(feature = "index-binned")]
pub struct BinnedIndex;

#[cfg(feature = "index-binned")]
impl IntervalIndex for BinnedIndex {
    fn intersect_pairs(
        a: &[Transcript],
        b: &[Transcript],
        opts: &IntersectOpts,
    ) -> Vec<(usize, usize)> {
        super::index_binned::binned_intersect_pairs(a, b, opts)
    }
}

/// Intersect pairs using the sweep-line backend.
///
/// Prefer this for already-sorted inputs; it is typically faster than indexed approaches when both
/// sides are sorted.
pub fn intersect_pairs(
    a: &[Transcript],
    b: &[Transcript],
    opts: &IntersectOpts,
) -> Vec<(usize, usize)> {
    SweepIndex::intersect_pairs(a, b, opts)
}

/// Intersect pairs using the binned index backend.
#[cfg(feature = "index-binned")]
pub fn intersect_pairs_binned(
    a: &[Transcript],
    b: &[Transcript],
    opts: &IntersectOpts,
) -> Vec<(usize, usize)> {
    BinnedIndex::intersect_pairs(a, b, opts)
}
