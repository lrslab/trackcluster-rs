use std::collections::HashMap;

use crate::model::{Interval, Transcript};

use super::{partition, IntersectOpts, PartitionKey, StrandMode};

const BIN_SIZE: u32 = 16_384;

fn span(transcript: &Transcript) -> Interval {
    Interval {
        start: transcript.tx_start,
        end: transcript.tx_end,
    }
}

fn span_overlap_len(a: &Transcript, b: &Transcript) -> u32 {
    span(a).overlap_len(span(b))
}

struct BinnedPartitionIndex {
    max_bin: u32,
    bins: Vec<Vec<usize>>,
}

impl BinnedPartitionIndex {
    fn build(b: &[Transcript], b_indices: &[usize]) -> Self {
        let max_end = b_indices
            .iter()
            .map(|&idx| b[idx].tx_end.get())
            .max()
            .unwrap_or(0);
        let max_bin = max_end / BIN_SIZE;

        let mut bins: Vec<Vec<usize>> = vec![Vec::new(); (max_bin + 1) as usize];
        for &idx in b_indices {
            let start = b[idx].tx_start.get();
            let end = b[idx].tx_end.get();
            if start >= end {
                continue;
            }

            let start_bin = start / BIN_SIZE;
            let end_bin = end.saturating_sub(1) / BIN_SIZE;
            for bin in start_bin..=end_bin {
                bins[bin as usize].push(idx);
            }
        }

        Self { max_bin, bins }
    }

    fn collect_candidates(
        &self,
        query_start: u32,
        query_end: u32,
        seen: &mut [u32],
        stamp: u32,
        out: &mut Vec<usize>,
    ) {
        if query_start >= query_end {
            return;
        }

        let start_bin = query_start / BIN_SIZE;
        if start_bin > self.max_bin {
            return;
        }

        let end_bin = (query_end.saturating_sub(1) / BIN_SIZE).min(self.max_bin);
        for bin in start_bin..=end_bin {
            for &idx in &self.bins[bin as usize] {
                if seen[idx] != stamp {
                    seen[idx] = stamp;
                    out.push(idx);
                }
            }
        }
    }
}

#[derive(Clone, Debug)]
/// Reusable allocation scratch space for binned intersection queries.
pub struct BinnedIntersectScratch {
    seen: Vec<u32>,
    stamp: u32,
    candidates: Vec<usize>,
}

impl BinnedIntersectScratch {
    /// Allocate scratch space for a catalog with `b_len` records.
    pub fn new(b_len: usize) -> Self {
        Self {
            seen: vec![0; b_len],
            stamp: 0,
            candidates: Vec::new(),
        }
    }

    fn ensure_b_len(&mut self, b_len: usize) {
        if self.seen.len() != b_len {
            self.seen.clear();
            self.seen.resize(b_len, 0);
            self.stamp = 0;
        }
        self.candidates.clear();
    }

    fn next_stamp(&mut self) -> u32 {
        self.stamp = self.stamp.wrapping_add(1);
        if self.stamp == 0 {
            self.seen.fill(0);
            self.stamp = 1;
        }
        self.stamp
    }
}

/// Reusable fixed-bin index bound to one transcript catalog.
///
/// Queries must use the same `b` records, in the same order, that were supplied to
/// [`BinnedIntersectIndex::build`]. The index stores catalog indices and partition/bin membership,
/// but does not retain or borrow `b`, so it cannot verify record identity at query time.
pub struct BinnedIntersectIndex {
    strand_mode: StrandMode,
    b_len: usize,
    partitions: HashMap<PartitionKey, BinnedPartitionIndex>,
}

impl BinnedIntersectIndex {
    /// Build an index for `b` using the selected strand policy.
    ///
    /// The returned index is valid only for queries using this exact catalog and record order.
    pub fn build(b: &[Transcript], strand_mode: StrandMode) -> Self {
        let b_parts = partition(b, strand_mode);
        let mut partitions: HashMap<PartitionKey, BinnedPartitionIndex> =
            HashMap::with_capacity(b_parts.len());
        for (key, indices) in b_parts {
            partitions.insert(key, BinnedPartitionIndex::build(b, &indices));
        }
        Self {
            strand_mode,
            b_len: b.len(),
            partitions,
        }
    }

    /// Replace the contents of `out` with matching pairs, reusing caller-provided scratch
    /// allocation.
    ///
    /// `b` must contain the same records in the same order as the catalog used to build this
    /// index, and `opts.strand_mode` must match the build-time strand mode. This method checks
    /// only the strand mode and catalog length. A detected mismatch clears `out` and returns no
    /// pairs; a same-length catalog or ordering mismatch is not detected and can produce
    /// incomplete or otherwise meaningless results.
    pub fn intersect_pairs_into(
        &self,
        a: &[Transcript],
        b: &[Transcript],
        opts: &IntersectOpts,
        scratch: &mut BinnedIntersectScratch,
        out: &mut Vec<(usize, usize)>,
    ) {
        out.clear();

        if opts.strand_mode != self.strand_mode || b.len() != self.b_len {
            return;
        }

        scratch.ensure_b_len(b.len());

        let a_parts = partition(a, opts.strand_mode);
        for (key, a_indices) in a_parts {
            let Some(index) = self.partitions.get(&key) else {
                continue;
            };

            for &a_idx in &a_indices {
                let a_tx = &a[a_idx];
                let start = a_tx.tx_start.get();
                let end = a_tx.tx_end.get();
                if start >= end {
                    continue;
                }

                let stamp = scratch.next_stamp();
                scratch.candidates.clear();
                index.collect_candidates(
                    start,
                    end,
                    &mut scratch.seen,
                    stamp,
                    &mut scratch.candidates,
                );

                for &b_idx in &scratch.candidates {
                    let b_tx = &b[b_idx];
                    let overlap = span_overlap_len(a_tx, b_tx);
                    if overlap == 0 {
                        continue;
                    }
                    if let Some(min_overlap) = opts.min_overlap_bp {
                        if overlap < min_overlap {
                            continue;
                        }
                    }
                    out.push((a_idx, b_idx));
                }
            }
        }

        out.sort_unstable();
    }

    /// Return matching pairs for a single query collection.
    ///
    /// The `b` catalog and strand-mode requirements are the same as for
    /// [`BinnedIntersectIndex::intersect_pairs_into`]. A detected strand-mode or length mismatch
    /// returns an empty vector.
    pub fn intersect_pairs(
        &self,
        a: &[Transcript],
        b: &[Transcript],
        opts: &IntersectOpts,
    ) -> Vec<(usize, usize)> {
        let mut scratch = BinnedIntersectScratch::new(b.len());
        let mut out = Vec::new();
        self.intersect_pairs_into(a, b, opts, &mut scratch, &mut out);
        out
    }
}

/// Returns all `(a_index, b_index)` pairs whose transcript spans overlap (half-open).
///
/// This backend builds a fixed-size binned index for each partition of `b` and queries it for each
/// transcript in `a`.
pub fn binned_intersect_pairs(
    a: &[Transcript],
    b: &[Transcript],
    opts: &IntersectOpts,
) -> Vec<(usize, usize)> {
    let index = BinnedIntersectIndex::build(b, opts.strand_mode);
    index.intersect_pairs(a, b, opts)
}

#[cfg(test)]
mod tests {
    use crate::interval::sort::sort_by_coord;
    use crate::interval::{sweep_intersect_pairs, StrandMode};
    use crate::model::{Bed12Attrs, Coord, Interval, Strand, Transcript};
    use proptest::prelude::*;

    use super::*;

    fn make_tx(chrom: String, strand: Strand, start: u32, len: u32, name: String) -> Transcript {
        let end = start.saturating_add(len);
        Transcript::new(
            chrom,
            strand,
            Coord::new(start),
            Coord::new(end),
            name,
            vec![Interval::new(Coord::new(start), Coord::new(end)).unwrap()],
            Bed12Attrs {
                score: 0,
                thick_start: Coord::new(start),
                thick_end: Coord::new(end),
                item_rgb: "0".to_owned(),
                extra_fields: Vec::new(),
            },
        )
        .unwrap()
    }

    fn naive_pairs(
        a: &[Transcript],
        b: &[Transcript],
        opts: &IntersectOpts,
    ) -> Vec<(usize, usize)> {
        let mut out = Vec::new();
        for (ai, atx) in a.iter().enumerate() {
            for (bi, btx) in b.iter().enumerate() {
                if atx.chrom != btx.chrom {
                    continue;
                }
                if opts.strand_mode == StrandMode::Match && atx.strand != btx.strand {
                    continue;
                }
                let overlap = span_overlap_len(atx, btx);
                if overlap == 0 {
                    continue;
                }
                if let Some(min_overlap) = opts.min_overlap_bp {
                    if overlap < min_overlap {
                        continue;
                    }
                }
                out.push((ai, bi));
            }
        }
        out.sort_unstable();
        out
    }

    proptest! {
        #[test]
        fn binned_matches_sweep_and_naive(
            a_specs in prop::collection::vec((prop_oneof![Just("chr1".to_owned()), Just("chr2".to_owned())],
                                              prop_oneof![Just(Strand::Plus), Just(Strand::Minus), Just(Strand::Unknown)],
                                              0u32..200_000,
                                              1u32..50_000), 0..12),
            b_specs in prop::collection::vec((prop_oneof![Just("chr1".to_owned()), Just("chr2".to_owned())],
                                              prop_oneof![Just(Strand::Plus), Just(Strand::Minus), Just(Strand::Unknown)],
                                              0u32..200_000,
                                              1u32..50_000), 0..12),
            strand_mode in prop_oneof![Just(StrandMode::Ignore), Just(StrandMode::Match)],
            min_overlap_bp in prop_oneof![Just(None), (1u32..25).prop_map(Some)],
        ) {
            let mut a: Vec<Transcript> = a_specs.into_iter().enumerate().map(|(i, (chrom, strand, start, len))| {
                make_tx(chrom, strand, start, len, format!("a{i}"))
            }).collect();
            let mut b: Vec<Transcript> = b_specs.into_iter().enumerate().map(|(i, (chrom, strand, start, len))| {
                make_tx(chrom, strand, start, len, format!("b{i}"))
            }).collect();

            sort_by_coord(&mut a);
            sort_by_coord(&mut b);

            let opts = IntersectOpts {
                strand_mode,
                min_overlap_bp,
            };

            let binned = binned_intersect_pairs(&a, &b, &opts);
            let sweep = sweep_intersect_pairs(&a, &b, &opts);
            let naive = naive_pairs(&a, &b, &opts);

            prop_assert_eq!(&binned, &sweep);
            prop_assert_eq!(&binned, &naive);
        }
    }
}
