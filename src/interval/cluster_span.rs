use crate::model::{Interval, Transcript};

use super::StrandMode;

#[derive(Clone, Debug, PartialEq, Eq)]
/// One connected component of overlapping transcript spans.
pub struct RangeCluster {
    /// Reference sequence shared by all members.
    pub chrom: String,
    /// Strand key when strand matching is enabled.
    pub strand: Option<crate::model::Strand>,
    /// Bounding span of every member.
    pub span: Interval,
    /// Indices into the original input slice.
    pub members: Vec<usize>,
}

fn span(transcript: &Transcript) -> Interval {
    Interval {
        start: transcript.tx_start,
        end: transcript.tx_end,
    }
}

/// Cluster transcript spans from arbitrarily ordered input.
///
/// Half-open spans that merely touch are placed in separate clusters.
pub fn cluster_by_span(records: &[Transcript], strand_mode: StrandMode) -> Vec<RangeCluster> {
    let mut ordered_indices: Vec<usize> = (0..records.len()).collect();
    ordered_indices.sort_by(|&left_index, &right_index| {
        let left = &records[left_index];
        let right = &records[right_index];
        left.chrom
            .cmp(&right.chrom)
            .then_with(|| {
                strand_mode
                    .key_strand(left.strand)
                    .cmp(&strand_mode.key_strand(right.strand))
            })
            .then_with(|| left.tx_start.cmp(&right.tx_start))
            .then_with(|| left.tx_end.cmp(&right.tx_end))
            .then_with(|| left_index.cmp(&right_index))
    });

    let mut clusters: Vec<RangeCluster> = Vec::new();

    let mut current: Option<RangeCluster> = None;
    for index in ordered_indices {
        let transcript = &records[index];
        let tx_span = span(transcript);
        let key_strand = strand_mode.key_strand(transcript.strand);

        match current.as_mut() {
            None => {
                current = Some(RangeCluster {
                    chrom: transcript.chrom.clone(),
                    strand: key_strand,
                    span: tx_span,
                    members: vec![index],
                });
            }
            Some(cluster) => {
                if transcript.chrom != cluster.chrom || key_strand != cluster.strand {
                    clusters.push(current.take().unwrap());
                    current = Some(RangeCluster {
                        chrom: transcript.chrom.clone(),
                        strand: key_strand,
                        span: tx_span,
                        members: vec![index],
                    });
                    continue;
                }

                if tx_span.start < cluster.span.end {
                    if tx_span.end > cluster.span.end {
                        cluster.span.end = tx_span.end;
                    }
                    cluster.members.push(index);
                } else {
                    clusters.push(current.take().unwrap());
                    current = Some(RangeCluster {
                        chrom: transcript.chrom.clone(),
                        strand: key_strand,
                        span: tx_span,
                        members: vec![index],
                    });
                }
            }
        }
    }

    if let Some(cluster) = current {
        clusters.push(cluster);
    }

    clusters.sort_by(|left, right| {
        left.chrom
            .cmp(&right.chrom)
            .then_with(|| left.span.start.cmp(&right.span.start))
            .then_with(|| left.span.end.cmp(&right.span.end))
            .then_with(|| left.strand.cmp(&right.strand))
            .then_with(|| left.members.cmp(&right.members))
    });

    clusters
}

#[cfg(test)]
mod tests {
    use crate::interval::sort::sort_by_coord;
    use crate::model::{Bed12Attrs, Coord, Interval, Strand, Transcript};

    use super::*;

    fn make_tx(chrom: &str, strand: Strand, start: u32, end: u32, name: &str) -> Transcript {
        Transcript::new(
            chrom.to_owned(),
            strand,
            Coord::new(start),
            Coord::new(end),
            name.to_owned(),
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

    #[test]
    fn cluster_by_span_groups_overlapping_records() {
        let mut records = vec![
            make_tx("chr1", Strand::Plus, 10, 20, "a"),
            make_tx("chr1", Strand::Plus, 15, 25, "b"),
            make_tx("chr1", Strand::Plus, 30, 40, "c"),
        ];
        sort_by_coord(&mut records);

        let clusters = cluster_by_span(&records, StrandMode::Ignore);
        assert_eq!(clusters.len(), 2);
        assert_eq!(clusters[0].members.len(), 2);
        assert_eq!(clusters[1].members.len(), 1);
    }

    #[test]
    fn cluster_by_span_respects_strand_mode() {
        let mut records = vec![
            make_tx("chr1", Strand::Plus, 10, 20, "a"),
            make_tx("chr1", Strand::Minus, 15, 25, "b"),
        ];
        sort_by_coord(&mut records);

        let clusters = cluster_by_span(&records, StrandMode::Ignore);
        assert_eq!(clusters.len(), 1);

        let clusters = cluster_by_span(&records, StrandMode::Match);
        assert_eq!(clusters.len(), 2);
    }

    #[test]
    fn cluster_by_span_handles_unsorted_strand_interleaving() {
        let records = vec![
            make_tx("chr1", Strand::Plus, 40, 60, "plus-right"),
            make_tx("chr1", Strand::Minus, 20, 30, "minus"),
            make_tx("chr1", Strand::Plus, 10, 50, "plus-left"),
        ];

        let clusters = cluster_by_span(&records, StrandMode::Match);
        assert_eq!(clusters.len(), 2);
        let plus = clusters
            .iter()
            .find(|cluster| cluster.strand == Some(Strand::Plus))
            .unwrap();
        assert_eq!(
            plus.span,
            Interval::new(Coord::new(10), Coord::new(60)).unwrap()
        );
        assert_eq!(plus.members, vec![2, 0]);
    }
}
