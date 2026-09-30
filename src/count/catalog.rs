//! Direct nearest-isoform assignment against an immutable catalog.

use std::collections::{BTreeMap, HashMap};
use std::io::Write;
use std::path::{Path, PathBuf};

use super::{
    assignment_score, gene_metadata_compatible, gene_name, next_catalog_stamp,
    read_structures_by_name, transcript_exon_overlap, AssignmentScore, SpanCandidateIndex,
    UniqueAssignmentOptions,
};
use crate::flow::artifact_manifest::atomic_write_with;
use crate::model::Transcript;

/// Why a read has no candidate in the fixed catalog.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum UnassignedReason {
    /// No catalog transcript overlaps the read's genomic span on its chromosome.
    NoCandidateLocus,
    /// Overlapping transcripts all have a different strand.
    StrandMismatch,
    /// Same-strand transcripts all have incompatible annotated gene IDs.
    GeneMismatch,
    /// Candidate spans overlap, but none of their exons overlap a read exon.
    NoExonOverlap,
}

impl UnassignedReason {
    fn as_str(self) -> &'static str {
        match self {
            Self::NoCandidateLocus => "no_candidate_locus",
            Self::StrandMismatch => "strand_mismatch",
            Self::GeneMismatch => "gene_mismatch",
            Self::NoExonOverlap => "no_exon_overlap",
        }
    }
}

/// A molecule that cannot be assigned within the supplied catalog.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct UnassignedRead {
    /// Original read ID, or `sample::read` in manifest mode.
    pub read_id: String,
    /// Read gene annotation, or `none` when absent.
    pub gene: String,
    /// First candidate-domain requirement that could not be satisfied.
    pub reason: UnassignedReason,
}

/// Complete accounting of the input molecules, sorted by read ID.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CatalogAssignmentResult {
    /// Number of input BED records, including identical duplicate alignments.
    pub input_records: usize,
    /// Exactly one pair per assigned molecule; no catalog records are created.
    pub read_to_isoform: Vec<(String, String)>,
    /// One diagnostic per unassigned molecule.
    pub unassigned_reads: Vec<UnassignedRead>,
}

/// Assign all input molecules to the nearest overlapping catalog isoform.
///
/// Candidates must share chromosome and strand, have compatible gene metadata
/// when both sides are annotated, and overlap at least one exonic base. Splice
/// differences are ranked, not rejected: every read with a candidate gets one
/// assignment. The existing unique-assignment structural score is minimized,
/// with isoform ID as the final tie-break. Catalog coverage, reference status,
/// embedded memberships and pre-existing mapping files never affect selection.
pub fn assign_reads_to_catalog(
    reads: &[Transcript],
    isoforms: &[Transcript],
    options: UniqueAssignmentOptions,
) -> anyhow::Result<CatalogAssignmentResult> {
    anyhow::ensure!(!isoforms.is_empty(), "fixed isoform catalog is empty");
    crate::identity::validate_isoform_ids(isoforms)?;
    let reads_by_name: BTreeMap<_, _> = read_structures_by_name(reads)?.into_iter().collect();

    let mut by_chrom: HashMap<&str, Vec<usize>> = HashMap::new();
    for (idx, isoform) in isoforms.iter().enumerate() {
        by_chrom.entry(&isoform.chrom).or_default().push(idx);
    }
    let indices: HashMap<_, _> = by_chrom
        .into_iter()
        .map(|(chrom, ids)| (chrom, SpanCandidateIndex::build(isoforms, &ids)))
        .collect();
    let mut seen = vec![0u32; isoforms.len()];
    let mut stamp = 0;
    let mut candidates = Vec::new();
    let mut result = CatalogAssignmentResult {
        input_records: reads.len(),
        read_to_isoform: Vec::new(),
        unassigned_reads: Vec::new(),
    };

    for (read_id, read) in reads_by_name {
        candidates.clear();
        if let Some(index) = indices.get(read.chrom.as_str()) {
            let current_stamp = next_catalog_stamp(&mut seen, &mut stamp);
            index.collect_overlapping(read, isoforms, &mut seen, current_stamp, &mut candidates);
        }
        let mut reason = UnassignedReason::NoCandidateLocus;
        if !candidates.is_empty() {
            reason = UnassignedReason::StrandMismatch;
            candidates.retain(|&idx| isoforms[idx].strand == read.strand);
        }
        if !candidates.is_empty() {
            reason = UnassignedReason::GeneMismatch;
            candidates.retain(|&idx| gene_metadata_compatible(read, &isoforms[idx]));
        }
        if !candidates.is_empty() {
            reason = UnassignedReason::NoExonOverlap;
            candidates.retain(|&idx| transcript_exon_overlap(read, &isoforms[idx]) > 0);
        }

        let best: Option<(AssignmentScore, &str)> = candidates
            .iter()
            .map(|&idx| {
                let isoform = &isoforms[idx];
                (
                    assignment_score(read, isoform, options.junction_offset),
                    isoform.name.as_str(),
                )
            })
            .min();
        if let Some((_, isoform_id)) = best {
            result
                .read_to_isoform
                .push((read_id.to_owned(), isoform_id.to_owned()));
        } else {
            result.unassigned_reads.push(UnassignedRead {
                read_id: read_id.to_owned(),
                gene: gene_name(read).unwrap_or("none").to_owned(),
                reason,
            });
        }
    }
    Ok(result)
}

/// Diagnostic files accompanying direct catalog counts.
#[derive(Clone, Debug)]
pub struct CatalogOutputPaths {
    /// Headerless read-to-isoform TSV, compatible with existing mapping readers.
    pub mapping: PathBuf,
    /// Unassigned molecules and reasons.
    pub unassigned: PathBuf,
    /// Record/molecule counts and assigned/unassigned totals.
    pub stats: PathBuf,
    /// Effective assignment policy and junction tolerance.
    pub provenance: PathBuf,
}

impl CatalogOutputPaths {
    /// Use the same stem as a single-sample count CSV.
    pub fn for_count(count_csv: &Path) -> Self {
        Self {
            mapping: count_csv.with_extension("read_to_isoform.tsv"),
            unassigned: count_csv.with_extension("unassigned_reads.tsv"),
            stats: count_csv.with_extension("assignment_stats.tsv"),
            provenance: count_csv.with_extension("provenance.tsv"),
        }
    }

    /// Append suffixes to a multi-sample output prefix, preserving any dots.
    pub fn for_count_multi(prefix: &Path) -> Self {
        let append = |suffix: &str| {
            let mut path = prefix.as_os_str().to_os_string();
            path.push(suffix);
            PathBuf::from(path)
        };
        Self {
            mapping: append(".read_to_isoform.tsv"),
            unassigned: append(".unassigned_reads.tsv"),
            stats: append(".assignment_stats.tsv"),
            provenance: append(".unique_assignment.provenance.tsv"),
        }
    }

    /// Paths for the CLI's input/output alias preflight.
    pub fn labeled_paths(&self) -> [(&'static str, &Path); 4] {
        [
            ("catalog assignment mapping output", &self.mapping),
            ("unassigned-read output", &self.unassigned),
            ("assignment-statistics output", &self.stats),
            ("assignment-provenance output", &self.provenance),
        ]
    }

    /// Atomically replace each report after assignment succeeds.
    pub fn write(
        &self,
        result: &CatalogAssignmentResult,
        options: UniqueAssignmentOptions,
    ) -> anyhow::Result<()> {
        atomic_write_with(&self.mapping, |writer| {
            crate::cluster::output::write_read_to_isoform_tsv_writer(
                writer,
                &result.read_to_isoform,
            )
            .map_err(Into::into)
        })?;
        atomic_write_with(&self.unassigned, |writer| {
            let mut tsv = csv::WriterBuilder::new()
                .delimiter(b'\t')
                .from_writer(writer);
            tsv.write_record(["read_id", "gene", "reason"])?;
            for read in &result.unassigned_reads {
                tsv.write_record([&read.read_id, &read.gene, read.reason.as_str()])?;
            }
            tsv.flush()?;
            Ok(())
        })?;
        atomic_write_with(&self.stats, |writer| {
            let assigned = result.read_to_isoform.len();
            let unassigned = result.unassigned_reads.len();
            writeln!(writer, "metric\tcount")?;
            writeln!(writer, "input_records\t{}", result.input_records)?;
            writeln!(writer, "total_reads\t{}", assigned + unassigned)?;
            writeln!(writer, "assigned_reads\t{assigned}")?;
            writeln!(writer, "unassigned_reads\t{unassigned}")?;
            Ok(())
        })?;
        atomic_write_with(&self.provenance, |writer| {
            super::write_unique_assignment_provenance_to_writer(writer, options)?;
            writeln!(writer, "assignment_source\tfixed_catalog")?;
            writeln!(writer, "read_scope\tall_input_reads")?;
            writeln!(
                writer,
                "candidate_policy\tsame_chrom_strand_gene_exon_overlap"
            )?;
            writeln!(writer, "splice_conflicts\trank_not_reject")?;
            writeln!(writer, "tie_break\tisoform_id")?;
            Ok(())
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{Bed12Attrs, Coord, Interval, Strand};

    fn tx(name: &str, exons: &[(u32, u32)]) -> Transcript {
        let start = Coord::new(exons[0].0);
        let end = Coord::new(exons.last().unwrap().1);
        Transcript::new(
            "chr1".to_owned(),
            Strand::Plus,
            start,
            end,
            name.to_owned(),
            exons
                .iter()
                .map(|&(s, e)| Interval::new(Coord::new(s), Coord::new(e)).unwrap())
                .collect(),
            Bed12Attrs {
                score: 0,
                thick_start: start,
                thick_end: end,
                item_rgb: "0".to_owned(),
                extra_fields: Vec::new(),
            },
        )
        .unwrap()
    }

    fn assign(reads: &[Transcript], isoforms: &[Transcript]) -> CatalogAssignmentResult {
        assign_reads_to_catalog(reads, isoforms, UniqueAssignmentOptions::default()).unwrap()
    }

    #[test]
    fn ranks_splicing_before_ends_and_assigns_all_reads_without_memberships() {
        let isoforms = vec![
            tx("a_wrong_splice", &[(100, 220), (300, 400)]),
            tx("z_exact", &[(100, 180), (300, 400)]),
            tx("longer", &[(100, 180), (300, 450)]),
        ];
        let result = assign(
            &[
                tx("full", &[(100, 180), (300, 400)]),
                tx("truncated", &[(130, 180), (300, 395)]),
                tx("novel_junction", &[(100, 250), (320, 400)]),
            ],
            &isoforms,
        );
        assert_eq!(
            result.read_to_isoform,
            vec![
                ("full".into(), "z_exact".into()),
                ("novel_junction".into(), "a_wrong_splice".into()),
                ("truncated".into(), "z_exact".into()),
            ]
        );
        assert!(result.unassigned_reads.is_empty());
    }

    #[test]
    fn junction_tolerance_changes_the_nearest_assignment() {
        let isoforms = vec![
            tx("near_ends", &[(100, 200), (300, 400)]),
            tx("exact_splice", &[(100, 205), (305, 440)]),
        ];
        let reads = vec![tx("read", &[(100, 205), (305, 400)])];
        let exact = assign_reads_to_catalog(
            &reads,
            &isoforms,
            UniqueAssignmentOptions { junction_offset: 0 },
        )
        .unwrap();
        assert_eq!(exact.read_to_isoform[0].1, "exact_splice");
        assert_eq!(assign(&reads, &isoforms).read_to_isoform[0].1, "near_ends");
    }

    #[test]
    fn single_exon_reads_can_match_internal_exons_but_not_intronic_only_spans() {
        let isoforms = vec![tx("multi_exon", &[(100, 150), (200, 250), (300, 350)])];
        let result = assign(
            &[tx("internal", &[(210, 230)]), tx("intronic", &[(160, 180)])],
            &isoforms,
        );
        assert_eq!(
            result.read_to_isoform,
            vec![("internal".into(), "multi_exon".into())]
        );
        assert_eq!(result.unassigned_reads[0].read_id, "intronic");
        assert_eq!(
            result.unassigned_reads[0].reason,
            UnassignedReason::NoExonOverlap
        );
    }

    #[test]
    fn reports_locus_strand_and_gene_mismatches_and_accepts_gene_sets() {
        let mut isoform = tx("catalog", &[(100, 200)]);
        isoform.metadata_mut().set_gene_id("GENEA");
        let mut foreign_chrom = tx("foreign_chrom", &[(100, 200)]);
        foreign_chrom.chrom = "chr2".into();
        let mut opposite = tx("opposite", &[(100, 200)]);
        opposite.strand = Strand::Minus;
        let mut unknown_strand = tx("unknown_strand", &[(100, 200)]);
        unknown_strand.strand = Strand::Unknown;
        let mut wrong_gene = tx("wrong_gene", &[(100, 200)]);
        wrong_gene.metadata_mut().set_gene_id("GENEB");
        let mut gene_set = tx("gene_set", &[(100, 200)]);
        gene_set.metadata_mut().set_gene_id("GENEB||GENEA");
        let result = assign(
            &[
                foreign_chrom,
                opposite,
                unknown_strand,
                wrong_gene,
                gene_set,
                tx("unannotated", &[(100, 200)]),
                tx("distant", &[(900, 1000)]),
            ],
            &[isoform],
        );
        assert_eq!(
            result.read_to_isoform,
            vec![
                ("gene_set".into(), "catalog".into()),
                ("unannotated".into(), "catalog".into()),
            ]
        );
        let reasons: Vec<_> = result
            .unassigned_reads
            .iter()
            .map(|r| (r.read_id.as_str(), r.reason))
            .collect();
        assert_eq!(
            reasons,
            vec![
                ("distant", UnassignedReason::NoCandidateLocus),
                ("foreign_chrom", UnassignedReason::NoCandidateLocus),
                ("opposite", UnassignedReason::StrandMismatch),
                ("unknown_strand", UnassignedReason::StrandMismatch),
                ("wrong_gene", UnassignedReason::GeneMismatch),
            ]
        );
    }

    #[test]
    fn annotated_reads_can_use_an_unannotated_catalog_and_choose_once_across_genes() {
        let mut read = tx("read", &[(100, 200)]);
        read.metadata_mut().set_gene_id("GENEA");
        assert_eq!(
            assign(&[read], &[tx("plain_bed", &[(100, 200)])]).read_to_isoform[0].1,
            "plain_bed"
        );

        let mut far = tx("gene_a", &[(100, 230)]);
        far.metadata_mut().set_gene_id("GENEA");
        let mut near = tx("gene_b", &[(100, 200)]);
        near.metadata_mut().set_gene_id("GENEB");
        assert_eq!(
            assign(&[tx("read", &[(100, 200)])], &[far, near]).read_to_isoform,
            vec![("read".into(), "gene_b".into())]
        );
    }

    #[test]
    fn minus_strand_terminal_variants_choose_the_nearest_ends() {
        let mut isoforms = vec![
            tx("long", &[(100, 200), (300, 400)]),
            tx("short", &[(150, 200), (300, 400)]),
        ];
        for isoform in &mut isoforms {
            isoform.strand = Strand::Minus;
        }
        let mut read = tx("read", &[(152, 200), (300, 400)]);
        read.strand = Strand::Minus;
        assert_eq!(assign(&[read], &isoforms).read_to_isoform[0].1, "short");
    }

    #[test]
    fn ties_ignore_catalog_order_reference_status_and_embedded_read_support() {
        let mut reference = tx("z_reference", &[(100, 200)]);
        reference.metadata_mut().set_name2("read,ghost,|9999");
        reference.metadata_mut().set_transcript_type("isoform_anno");
        let novel = tx("a_novel", &[(100, 200)]);
        let reads = vec![tx("read", &[(100, 200)]), tx("read", &[(100, 200)])];
        let result = assign(&reads, &[reference.clone(), novel.clone()]);
        assert_eq!(result, assign(&reads, &[novel, reference]));
        assert_eq!(result.input_records, 2);
        assert_eq!(
            result.read_to_isoform,
            vec![("read".into(), "a_novel".into())]
        );
    }

    #[test]
    fn rejects_invalid_identity_inputs_and_accepts_empty_samples() {
        let catalog = vec![tx("iso", &[(100, 200)])];
        let duplicate_catalog = vec![catalog[0].clone(), catalog[0].clone()];
        assert!(assign_reads_to_catalog(
            &[],
            &duplicate_catalog,
            UniqueAssignmentOptions::default()
        )
        .unwrap_err()
        .to_string()
        .contains("duplicate isoform id"));
        let conflicting = vec![tx("read", &[(100, 200)]), tx("read", &[(110, 200)])];
        assert!(assign_reads_to_catalog(
            &conflicting,
            &catalog,
            UniqueAssignmentOptions::default()
        )
        .unwrap_err()
        .to_string()
        .contains("conflicting alignments"));
        assert!(
            assign_reads_to_catalog(&[], &[], UniqueAssignmentOptions::default())
                .unwrap_err()
                .to_string()
                .contains("catalog is empty")
        );
        let empty = assign(&[], &catalog);
        assert_eq!(empty.input_records, 0);
        assert!(empty.read_to_isoform.is_empty());
        assert!(empty.unassigned_reads.is_empty());
    }
}
