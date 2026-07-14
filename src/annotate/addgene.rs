use std::collections::{BTreeSet, HashMap};

use crate::interval::{intersect_pairs, sort_by_coord, IntersectOpts, StrandMode};
use crate::model::{Interval, Transcript};

#[derive(Clone, Copy, Debug)]
pub struct AddGeneOpts {
    pub fraction_read: f64,
    pub fraction_ref: f64,
}

impl Default for AddGeneOpts {
    fn default() -> Self {
        Self {
            fraction_read: 0.01,
            fraction_ref: 0.05,
        }
    }
}

impl AddGeneOpts {
    /// Validate overlap fractions before gene assignment.
    pub fn validate(self) -> Result<(), crate::config::ParameterError> {
        crate::config::UnitFraction::new("read overlap fraction", self.fraction_read)?;
        crate::config::UnitFraction::new("reference overlap fraction", self.fraction_ref)?;
        Ok(())
    }
}

fn exon_len(tx: &Transcript) -> u32 {
    tx.exons.iter().map(|exon| exon.len()).sum()
}

fn span(tx: &Transcript) -> Interval {
    Interval {
        start: tx.tx_start,
        end: tx.tx_end,
    }
}

fn span_len(tx: &Transcript) -> u32 {
    tx.tx_end.get().saturating_sub(tx.tx_start.get())
}

fn span_overlap_len(a: &Transcript, b: &Transcript) -> u32 {
    span(a).overlap_len(span(b))
}

fn gene_name(tx: &Transcript) -> &str {
    tx.metadata().gene_id_field().unwrap_or(tx.name.as_str())
}

pub fn dedup_longest_by_name(reads: &[Transcript]) -> Vec<Transcript> {
    let mut out: Vec<Transcript> = Vec::new();
    let mut pos: HashMap<String, usize> = HashMap::new();

    for read in reads {
        match pos.get(&read.name).copied() {
            None => {
                out.push(read.clone());
                pos.insert(read.name.clone(), out.len() - 1);
            }
            Some(existing_idx) => {
                let existing_len = exon_len(&out[existing_idx]);
                let candidate_len = exon_len(read);
                if candidate_len > existing_len {
                    out[existing_idx] = read.clone();
                }
            }
        }
    }

    out
}

fn add_gene_no_dedup_impl(
    reads: &[Transcript],
    references: &[Transcript],
    opts: AddGeneOpts,
) -> Vec<Transcript> {
    let mut reads = reads.to_vec();
    for read in &mut reads {
        read.metadata_mut().set_gene_id("none");
    }

    let mut reads_sorted = reads.clone();
    let mut refs_sorted: Vec<Transcript> = references.to_vec();
    sort_by_coord(&mut reads_sorted);
    sort_by_coord(&mut refs_sorted);

    let candidates = intersect_pairs(
        &reads_sorted,
        &refs_sorted,
        &IntersectOpts {
            strand_mode: StrandMode::Match,
            min_overlap_bp: None,
        },
    );

    let mut genes_by_read: HashMap<&str, BTreeSet<&str>> = HashMap::new();
    for (read_idx, ref_idx) in candidates {
        let read = &reads_sorted[read_idx];
        let reference = &refs_sorted[ref_idx];

        let overlap = span_overlap_len(read, reference);
        if overlap == 0 {
            continue;
        }

        let read_len = span_len(read);
        let ref_len = span_len(reference);
        if read_len == 0 || ref_len == 0 {
            continue;
        }

        let overlap_f = overlap as f64;
        if overlap_f / (read_len as f64) < opts.fraction_read {
            continue;
        }
        if overlap_f / (ref_len as f64) < opts.fraction_ref {
            continue;
        }

        genes_by_read
            .entry(read.name.as_str())
            .or_default()
            .insert(gene_name(reference));
    }

    for read in &mut reads {
        let Some(genes) = genes_by_read.get(read.name.as_str()) else {
            continue;
        };
        if genes.is_empty() {
            continue;
        }

        let joined = genes.iter().copied().collect::<Vec<_>>().join("||");
        read.metadata_mut().set_gene_id(joined);
    }

    reads
}

/// Assign genes without deduplicating reads, returning invalid option errors.
pub fn try_add_gene_no_dedup(
    reads: &[Transcript],
    references: &[Transcript],
    opts: AddGeneOpts,
) -> Result<Vec<Transcript>, crate::config::ParameterError> {
    opts.validate()?;
    Ok(add_gene_no_dedup_impl(reads, references, opts))
}

pub fn add_gene_no_dedup(
    reads: &[Transcript],
    references: &[Transcript],
    opts: AddGeneOpts,
) -> Vec<Transcript> {
    try_add_gene_no_dedup(reads, references, opts)
        .unwrap_or_else(|error| panic!("invalid gene-assignment options: {error}"))
}

/// Deduplicate and assign genes, returning invalid option errors.
pub fn try_add_gene(
    reads: &[Transcript],
    references: &[Transcript],
    opts: AddGeneOpts,
) -> Result<Vec<Transcript>, crate::config::ParameterError> {
    opts.validate()?;
    let reads = dedup_longest_by_name(reads);
    Ok(add_gene_no_dedup_impl(&reads, references, opts))
}

pub fn add_gene(
    reads: &[Transcript],
    references: &[Transcript],
    opts: AddGeneOpts,
) -> Vec<Transcript> {
    try_add_gene(reads, references, opts)
        .unwrap_or_else(|error| panic!("invalid gene-assignment options: {error}"))
}

#[cfg(test)]
mod tests {
    use crate::model::{Bed12Attrs, Coord, Interval, Strand, Transcript};

    use super::*;

    fn make_tx(
        name: &str,
        start: u32,
        end: u32,
        exons: &[(u32, u32)],
        gene_name: &str,
    ) -> Transcript {
        let exons = exons
            .iter()
            .map(|(s, e)| Interval::new(Coord::new(*s), Coord::new(*e)).unwrap())
            .collect::<Vec<_>>();

        Transcript::new(
            "chr1".to_owned(),
            Strand::Plus,
            Coord::new(start),
            Coord::new(end),
            name.to_owned(),
            exons,
            Bed12Attrs {
                score: 0,
                thick_start: Coord::new(start),
                thick_end: Coord::new(end),
                item_rgb: "0".to_owned(),
                extra_fields: vec![
                    "none".to_owned(),
                    "none".to_owned(),
                    "none".to_owned(),
                    "-1,".to_owned(),
                    "isoform_anno".to_owned(),
                    gene_name.to_owned(),
                    "none".to_owned(),
                    "none".to_owned(),
                ],
            },
        )
        .unwrap()
    }

    #[test]
    fn assigns_gene_name_when_overlap_passes_thresholds() {
        let refs = vec![make_tx("ref", 100, 200, &[(100, 150), (160, 200)], "GENE1")];
        let reads = vec![make_tx("read1", 120, 180, &[(120, 180)], "none")];

        let out = add_gene(&reads, &refs, AddGeneOpts::default());
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].metadata().gene_id(), Some("GENE1"));
    }

    #[test]
    fn deduplicates_reads_by_longest_exon_sum() {
        let refs = vec![make_tx("ref", 0, 1000, &[(0, 1000)], "GENE1")];
        let reads = vec![
            make_tx("dup", 100, 150, &[(100, 120), (130, 150)], "none"),
            make_tx("dup", 100, 200, &[(100, 200)], "none"),
        ];

        let out = add_gene(&reads, &refs, AddGeneOpts::default());
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].tx_end.get(), 200);
    }

    #[test]
    fn add_gene_no_dedup_preserves_duplicate_reads() {
        let refs = vec![make_tx("ref", 0, 1000, &[(0, 1000)], "GENE1")];
        let reads = vec![
            make_tx("dup", 100, 150, &[(100, 150)], "none"),
            make_tx("dup", 100, 200, &[(100, 200)], "none"),
        ];

        let out = add_gene_no_dedup(&reads, &refs, AddGeneOpts::default());
        assert_eq!(out.len(), 2);
        assert_eq!(out[0].metadata().gene_id(), Some("GENE1"));
        assert_eq!(out[1].metadata().gene_id(), Some("GENE1"));
    }

    #[test]
    fn add_gene_marks_unassigned_reads_as_none() {
        let reads = vec![make_tx("unmatched", 100, 200, &[(100, 200)], "STALE_GENE")];

        let out = add_gene_no_dedup(&reads, &[], AddGeneOpts::default());
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].metadata().gene_id_field(), Some("none"));
    }

    #[test]
    fn rejects_invalid_overlap_fractions_at_library_boundary() {
        for opts in [
            AddGeneOpts {
                fraction_read: f64::NAN,
                fraction_ref: 0.05,
            },
            AddGeneOpts {
                fraction_read: -0.1,
                fraction_ref: 0.05,
            },
            AddGeneOpts {
                fraction_read: 0.01,
                fraction_ref: 1.1,
            },
        ] {
            assert!(try_add_gene(&[], &[], opts).is_err());
            assert!(try_add_gene_no_dedup(&[], &[], opts).is_err());
        }
    }
}
