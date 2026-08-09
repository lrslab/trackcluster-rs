//! Strict Dorado/modBAM import into normalized modification observations.
//!
//! V1 emits explicit `MM`/`ML` calls and projects unlisted canonical bases as
//! low-probability or unknown observations according to the MM skip marker.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io;
use std::path::{Path, PathBuf};

use noodles_bam as bam;
use noodles_sam as sam;
use sam::alignment::record::cigar::{op::Kind, Op};
use sam::alignment::record::data::field::{value::Array, Tag, Value};
use sam::record::data::field::value::{
    base_modifications::{
        group::{Modification, Status, Strand as MmStrand, UnmodifiedBase},
        Group,
    },
    BaseModifications,
};
use thiserror::Error;

use crate::model::Strand;
use crate::modification::{
    ImplicitSkipPolicy, ModObservation, ModObservationKey, ModSiteKey, ObservationState,
};

/// Policy for a structurally invalid alignment record.
///
/// File-open and BAM framing/header errors always fail because continuing from
/// them is not guaranteed to be safe.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum InvalidRecordPolicy {
    /// Stop at the first invalid record.
    #[default]
    Fail,
    /// Exclude invalid records and count a stable reason in [`ModBamQc`].
    Skip,
}

/// Conversion used for SAM `ML:B:C` probability bytes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MlProbabilitySemantics {
    /// `N` denotes `[N/256, (N+1)/256)` (upper bound closed for 255),
    /// and the normalized point estimate is the interval midpoint.
    SamIntervalMidpoint,
}

impl MlProbabilitySemantics {
    /// Stable provenance token for assay metadata or audit output.
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::SamIntervalMidpoint => "sam_ml_interval_midpoint_v1",
        }
    }
}

/// Interpretation override for candidates omitted from an `MM` group marked `?`.
///
/// The SAM specification defines `?` as unknown. Some producers can also use
/// `?` while sparsifying a known all-context candidate universe by an emission
/// threshold. Reinterpreting those omissions therefore requires an explicit
/// source-specific assertion.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum MmQuestionMarkPolicy {
    /// Preserve the SAM meaning: omitted candidates have unknown state.
    #[default]
    Unknown,
    /// Treat omissions as known to be below a declared positive source threshold.
    BelowEmissionThreshold,
}

impl MmQuestionMarkPolicy {
    /// Stable provenance token for QC and assay metadata.
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Unknown => "unknown",
            Self::BelowEmissionThreshold => "below_emission_threshold",
        }
    }
}

/// Return the probability interval represented by a SAM `ML` byte.
pub const fn ml_probability_interval(value: u8) -> (f64, f64) {
    (value as f64 / 256.0, (value as f64 + 1.0) / 256.0)
}

/// Convert a SAM `ML` byte to the fixed V1 interval-midpoint estimate.
pub const fn ml_byte_to_probability(value: u8) -> f64 {
    (value as f64 + 0.5) / 256.0
}

/// Dorado/modBAM import settings.
#[derive(Clone, Debug, PartialEq)]
pub struct ModBamOptions {
    /// Assay compatibility stratum copied into every observation.
    pub assay_id: String,
    /// Biological sample name; raw BAM query names become `<sample>::<name>`.
    pub sample: String,
    /// One SAM modification code, such as `A+a`, `C+m`, or `A+17596`.
    pub target_mod_code: String,
    /// Candidate universe expanded by this importer.
    pub candidate_rule: String,
    /// Minimum retained MAPQ. Missing MAPQ is interpreted as 0.
    pub min_mapq: u8,
    /// Handling of malformed tags, projections, and duplicate primary reads.
    pub invalid_record_policy: InvalidRecordPolicy,
    /// Dorado `--modified-bases-threshold`, when known.
    pub source_emission_threshold: Option<f64>,
    /// Explicit interpretation of unlisted candidates in `MM` groups marked `?`.
    pub mm_question_mark_policy: MmQuestionMarkPolicy,
}

impl ModBamOptions {
    /// Create options for one sample, assay, and target modification.
    pub fn new(
        assay_id: impl Into<String>,
        sample: impl Into<String>,
        target_mod_code: impl Into<String>,
    ) -> Self {
        Self {
            assay_id: assay_id.into(),
            sample: sample.into(),
            target_mod_code: target_mod_code.into(),
            candidate_rule: "all-target-canonical-bases".to_owned(),
            min_mapq: 0,
            invalid_record_policy: InvalidRecordPolicy::Fail,
            source_emission_threshold: None,
            mm_question_mark_policy: MmQuestionMarkPolicy::Unknown,
        }
    }
}

impl Default for ModBamOptions {
    fn default() -> Self {
        Self::new("dorado", "sample", "A+a")
    }
}

/// Stable reason for excluding an invalid record under skip policy.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum InvalidRecordReason {
    /// The BAM record itself could not be decoded.
    RecordDecode,
    /// A mapped primary record has no usable query name.
    MissingReadName,
    /// A query name contains a control character.
    InvalidReadName,
    /// A mapped record has no valid reference dictionary entry.
    InvalidReference,
    /// A mapped record has no valid alignment start.
    InvalidAlignmentStart,
    /// An auxiliary field could not be decoded.
    InvalidAuxData,
    /// `MM` is absent.
    MissingMm,
    /// `MM` is duplicated.
    DuplicateMm,
    /// `MM` is not a `Z` string.
    InvalidMmType,
    /// `ML` is absent.
    MissingMl,
    /// `ML` is duplicated.
    DuplicateMl,
    /// `ML` is not a `B:C` unsigned-byte array.
    InvalidMlType,
    /// An `ML` array value could not be decoded.
    InvalidMlValue,
    /// `MN` is absent.
    MissingMn,
    /// `MN` is duplicated.
    DuplicateMn,
    /// `MN` is not a nonnegative integer representable as `usize`.
    InvalidMn,
    /// `MN` differs from the current decoded `SEQ` length.
    MnMismatch,
    /// noodles rejected the `MM` grammar or delta positions.
    InvalidMm,
    /// `ML` length does not equal all `MM` positions times all group codes.
    MlLengthMismatch,
    /// The CIGAR could not be decoded.
    InvalidCigar,
    /// CIGAR query/reference consumption is inconsistent or overflows.
    InvalidProjection,
    /// The CIGAR extends beyond its header reference length.
    ReferenceOutOfBounds,
    /// The target appears more than once at the same read/genomic site.
    DuplicateTargetObservation,
    /// More than one `MM` group describes the selected target code.
    DuplicateTargetGroup,
    /// An explicit target call falls outside the declared motif candidate universe.
    CandidateRuleMismatch,
    /// A generated normalized observation violated the shared schema.
    InvalidObservation,
    /// More than one mapped primary alignment has the same query name.
    DuplicatePrimary,
}

/// Query-base relationship to the genomic reference after walking CIGAR.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum QueryBaseProjection {
    /// The query base aligns to this zero-based reference base.
    Reference(u32),
    /// The query base belongs to an insertion and has no reference base.
    Insertion,
    /// The query base is soft clipped and has no reference base.
    SoftClip,
}

/// Complete projection of stored `SEQ` indexes through one CIGAR.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct QueryReferenceProjection {
    /// One projection entry for every decoded base in stored BAM `SEQ` order.
    pub query_bases: Vec<QueryBaseProjection>,
    /// Zero-based exclusive reference end after all reference-consuming ops.
    pub reference_end0: u32,
}

/// Error returned by the pure CIGAR projector.
#[derive(Clone, Debug, Error, PartialEq, Eq)]
pub enum ProjectionError {
    /// Zero-length CIGAR operations are rejected.
    #[error("CIGAR operation {operation_index} has zero length")]
    ZeroLength {
        /// Zero-based CIGAR operation index.
        operation_index: usize,
    },
    /// A query-consuming operation would exceed `SEQ`.
    #[error(
        "CIGAR query length exceeds SEQ at operation {operation_index}: {attempted} > {sequence_len}"
    )]
    QueryOverrun {
        /// Zero-based CIGAR operation index.
        operation_index: usize,
        /// Attempted zero-based exclusive query end.
        attempted: usize,
        /// Decoded BAM `SEQ` length.
        sequence_len: usize,
    },
    /// Total CIGAR query consumption differs from `SEQ` length.
    #[error("CIGAR consumes {cigar_query_len} query bases but SEQ has {sequence_len}")]
    QueryLengthMismatch {
        /// Total number of query bases consumed by CIGAR.
        cigar_query_len: usize,
        /// Decoded BAM `SEQ` length.
        sequence_len: usize,
    },
    /// Reference-coordinate arithmetic exceeded `u32`.
    #[error("CIGAR reference coordinate overflows u32 at operation {operation_index}")]
    ReferenceOverflow {
        /// Zero-based CIGAR operation index.
        operation_index: usize,
    },
    /// A mapped alignment contains no query base aligned to a reference base.
    #[error("mapped CIGAR contains no M, =, or X query base")]
    NoAlignedQueryBase,
}

/// Project decoded BAM `SEQ` indexes to zero-based reference bases.
///
/// `alignment_start0` is the zero-based alignment start. `M`, `=`, and `X`
/// consume both axes; `I` and `S` consume query only; `D` and `N` consume
/// reference only; `H` and `P` consume neither.
pub fn project_query_to_reference(
    alignment_start0: u32,
    sequence_len: usize,
    cigar: &[Op],
) -> Result<QueryReferenceProjection, ProjectionError> {
    let mut query_bases = Vec::with_capacity(sequence_len);
    let mut reference_pos0 = alignment_start0;
    let mut has_aligned_query_base = false;

    for (operation_index, operation) in cigar.iter().copied().enumerate() {
        let len = operation.len();
        if len == 0 {
            return Err(ProjectionError::ZeroLength { operation_index });
        }

        let checked_query_end = |current: usize| {
            current
                .checked_add(len)
                .filter(|&end| end <= sequence_len)
                .ok_or(ProjectionError::QueryOverrun {
                    operation_index,
                    attempted: current.saturating_add(len),
                    sequence_len,
                })
        };
        let advance_reference = |current: u32| {
            let len = u32::try_from(len)
                .map_err(|_| ProjectionError::ReferenceOverflow { operation_index })?;
            current
                .checked_add(len)
                .ok_or(ProjectionError::ReferenceOverflow { operation_index })
        };

        match operation.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                let query_end = checked_query_end(query_bases.len())?;
                let reference_end0 = advance_reference(reference_pos0)?;
                for pos0 in reference_pos0..reference_end0 {
                    query_bases.push(QueryBaseProjection::Reference(pos0));
                }
                debug_assert_eq!(query_bases.len(), query_end);
                reference_pos0 = reference_end0;
                has_aligned_query_base = true;
            }
            Kind::Insertion => {
                let query_end = checked_query_end(query_bases.len())?;
                query_bases.resize(query_end, QueryBaseProjection::Insertion);
            }
            Kind::SoftClip => {
                let query_end = checked_query_end(query_bases.len())?;
                query_bases.resize(query_end, QueryBaseProjection::SoftClip);
            }
            Kind::Deletion | Kind::Skip => {
                reference_pos0 = advance_reference(reference_pos0)?;
            }
            Kind::HardClip | Kind::Pad => {}
        }
    }

    if query_bases.len() != sequence_len {
        return Err(ProjectionError::QueryLengthMismatch {
            cigar_query_len: query_bases.len(),
            sequence_len,
        });
    }
    if !has_aligned_query_base {
        return Err(ProjectionError::NoAlignedQueryBase);
    }

    Ok(QueryReferenceProjection {
        query_bases,
        reference_end0: reference_pos0,
    })
}

/// Deterministic import and filtering counters.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ModBamQc {
    /// All BAM records encountered, including records later filtered.
    pub total_records: usize,
    /// Unmapped records excluded unconditionally.
    pub skipped_unmapped: usize,
    /// Secondary records excluded unconditionally.
    pub skipped_secondary: usize,
    /// Supplementary records excluded unconditionally.
    pub skipped_supplementary: usize,
    /// Unique primary records excluded by MAPQ.
    pub skipped_below_mapq: usize,
    /// Query names with two or more retained-class primary alignments.
    pub duplicate_primary_reads: usize,
    /// All primary records discarded because their query name was duplicated.
    pub skipped_duplicate_primary_records: usize,
    /// Invalid unique records excluded under skip policy, by stable reason.
    pub invalid_records: BTreeMap<InvalidRecordReason, usize>,
    /// Valid unique primary records whose tags and CIGAR were imported.
    pub retained_records: usize,
    /// Target canonical bases in retained records before candidate-rule filtering.
    pub target_canonical_bases: usize,
    /// Target bases satisfying the declared candidate rule.
    pub target_candidate_bases: usize,
    /// Retained records with no `MM` group containing the target code.
    pub records_without_target_group: usize,
    /// Number of target-containing groups with no explicit skip marker.
    pub target_groups_skip_flag_omitted: usize,
    /// Number of target-containing groups marked `.`.
    pub target_groups_low_probability: usize,
    /// Number of target-containing groups marked `?`.
    pub target_groups_unknown: usize,
    /// Unlisted target canonical bases interpreted as low probability.
    pub implicit_low_probability_candidates: usize,
    /// Unlisted target canonical bases whose state is unknown.
    pub unknown_candidates: usize,
    /// All `ML` bytes consumed, including non-target groups/codes.
    pub ml_values_consumed: usize,
    /// Explicit target calls decoded before genomic projection.
    pub explicit_target_calls: usize,
    /// Explicit target calls on inserted query bases.
    pub target_calls_in_insertions: usize,
    /// Explicit target calls on soft-clipped query bases.
    pub target_calls_in_soft_clips: usize,
    /// Implicit or unknown target candidates on inserted query bases.
    pub implicit_calls_in_insertions: usize,
    /// Implicit or unknown target candidates on soft-clipped query bases.
    pub implicit_calls_in_soft_clips: usize,
    /// Explicit target observations emitted with genomic coordinates.
    pub emitted_explicit_observations: usize,
    /// Low-probability implicit observations emitted with genomic coordinates.
    pub emitted_implicit_observations: usize,
    /// Unknown observations emitted with genomic coordinates.
    pub emitted_unknown_observations: usize,
    /// All target observations emitted with genomic coordinates.
    pub emitted_observations: usize,
}

/// Dataset-level semantics that must accompany imported observations.
#[derive(Clone, Debug, PartialEq)]
pub struct ModBamImportSemantics {
    /// Canonicalized single target code used for filtering.
    pub target_mod_code: String,
    /// Canonicalized all-context or centered-IUPAC candidate rule.
    pub candidate_rule: String,
    /// Declared Dorado emission threshold, if known.
    pub source_emission_threshold: Option<f64>,
    /// Applied interpretation of candidates omitted from `MM` groups marked `?`.
    pub mm_question_mark_policy: MmQuestionMarkPolicy,
    /// Combined interpretation of omitted candidates in target groups.
    pub implicit_skip_policy: ImplicitSkipPolicy,
    /// Whether the retained candidate universe is demonstrably complete.
    pub candidate_observations_complete: bool,
    /// Whether every emitted observation has an explicit ML probability.
    pub explicit_observations_only: bool,
    /// Exact conversion used for all normalized probabilities.
    pub ml_probability_semantics: MlProbabilitySemantics,
}

/// Imported normalized calls, QC counters, and interpretation provenance.
#[derive(Clone, Debug, PartialEq)]
pub struct ModBamImportResult {
    /// Explicit, genomically projected observations in deterministic key order.
    pub observations: Vec<ModObservation>,
    /// Filtering, validity, skip-semantic, and projection counts.
    pub qc: ModBamQc,
    /// Probability and candidate-completeness interpretation.
    pub semantics: ModBamImportSemantics,
}

/// Fatal modBAM import error.
#[derive(Debug, Error)]
pub enum ModBamError {
    /// Options violate the importer contract.
    #[error("invalid modBAM options: {0}")]
    InvalidOptions(String),
    /// Input could not be opened.
    #[error("open modBAM {path:?}")]
    Open {
        /// Input path.
        path: PathBuf,
        /// Underlying open error.
        #[source]
        source: io::Error,
    },
    /// BAM header framing or content is invalid.
    #[error("read modBAM header from {path:?}")]
    ReadHeader {
        /// Input path.
        path: PathBuf,
        /// Underlying BAM header error.
        #[source]
        source: io::Error,
    },
    /// A record is invalid and fail policy is active.
    #[error("invalid modBAM record {record_ordinal} ({read_name}): {reason:?}: {detail}")]
    InvalidRecord {
        /// One-based physical BAM record ordinal.
        record_ordinal: usize,
        /// Decoded query name, or `<unknown>` when unavailable.
        read_name: String,
        /// Stable QC category.
        reason: InvalidRecordReason,
        /// Record-specific diagnostic detail.
        detail: String,
    },
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct TargetModCode {
    unmodified_base: UnmodifiedBase,
    strand: MmStrand,
    modification: Modification,
    normalized: String,
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum CandidateRule {
    AllTargetCanonicalBases,
    CenteredIupacMotif {
        normalized: String,
        masks: Vec<u8>,
        center: usize,
    },
}

impl CandidateRule {
    fn parse(value: &str, target: &TargetModCode) -> Result<Self, String> {
        if value.eq_ignore_ascii_case("all-target-canonical-bases") {
            return Ok(Self::AllTargetCanonicalBases);
        }
        let normalized = crate::modification::types::normalize_centered_iupac_motif(
            value,
            u8::from(target.unmodified_base),
        )
        .map_err(|error| {
            error.replace(
                "must be an odd-length",
                "must be all-target-canonical-bases or an odd-length",
            )
        })?;
        let masks = normalized
            .bytes()
            .map(|base| {
                crate::modification::types::iupac_dna_mask(base).expect("normalized IUPAC motif")
            })
            .collect::<Vec<_>>();
        let center = masks.len() / 2;
        Ok(Self::CenteredIupacMotif {
            normalized,
            masks,
            center,
        })
    }

    fn normalized(&self) -> &str {
        match self {
            Self::AllTargetCanonicalBases => "all-target-canonical-bases",
            Self::CenteredIupacMotif { normalized, .. } => normalized,
        }
    }

    fn observation_context(&self) -> Option<&str> {
        match self {
            Self::AllTargetCanonicalBases => None,
            Self::CenteredIupacMotif { normalized, .. } => Some(normalized),
        }
    }

    fn candidate_positions(
        &self,
        decoded_sequence: &[u8],
        is_reverse_complemented: bool,
        target: &TargetModCode,
    ) -> Vec<usize> {
        match self {
            Self::AllTargetCanonicalBases => {
                if target.unmodified_base == UnmodifiedBase::N {
                    return (0..decoded_sequence.len()).collect();
                }
                let target_sequence_base = if is_reverse_complemented {
                    u8::from(target.unmodified_base.complement())
                } else {
                    u8::from(target.unmodified_base)
                };
                decoded_sequence
                    .iter()
                    .enumerate()
                    .filter_map(|(query_index, &base)| {
                        (base == target_sequence_base).then_some(query_index)
                    })
                    .collect()
            }
            Self::CenteredIupacMotif { masks, center, .. } => {
                let oriented_sequence = if is_reverse_complemented {
                    decoded_sequence
                        .iter()
                        .rev()
                        .map(|base| complement_sequence_base(*base))
                        .collect::<Vec<_>>()
                } else {
                    decoded_sequence
                        .iter()
                        .map(|base| base.to_ascii_uppercase())
                        .collect::<Vec<_>>()
                };
                let mut positions = oriented_sequence
                    .windows(masks.len())
                    .enumerate()
                    .filter(|(_, window)| {
                        window.iter().zip(masks).all(|(&base, &motif_mask)| {
                            sequence_base_mask(base).is_some_and(|mask| mask & motif_mask != 0)
                        })
                    })
                    .map(|(start, _)| {
                        let oriented_index = start + center;
                        if is_reverse_complemented {
                            decoded_sequence.len() - 1 - oriented_index
                        } else {
                            oriented_index
                        }
                    })
                    .collect::<Vec<_>>();
                positions.sort_unstable();
                positions
            }
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct TargetModel {
    target: TargetModCode,
    candidate_rule: CandidateRule,
}

fn sequence_base_mask(base: u8) -> Option<u8> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0b0001),
        b'C' => Some(0b0010),
        b'G' => Some(0b0100),
        b'T' | b'U' => Some(0b1000),
        _ => None,
    }
}

fn complement_sequence_base(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' | b'U' => b'A',
        b'R' => b'Y',
        b'Y' => b'R',
        b'S' => b'S',
        b'W' => b'W',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'D' => b'H',
        b'H' => b'D',
        b'V' => b'B',
        _ => b'N',
    }
}

fn parse_target_mod_code(value: &str) -> Result<TargetModCode, String> {
    let bytes = value.as_bytes();
    if bytes.len() < 3 {
        return Err(format!(
            "target_mod_code {value:?} must be a single code such as A+a or C+27551"
        ));
    }

    let unmodified_base = UnmodifiedBase::try_from(bytes[0])
        .map_err(|_| format!("target_mod_code {value:?} has an invalid canonical base"))?;
    let strand = MmStrand::try_from(bytes[1])
        .map_err(|_| format!("target_mod_code {value:?} has an invalid strand"))?;
    let raw_modification = &bytes[2..];
    let modification = if raw_modification.len() == 1 && raw_modification[0].is_ascii_lowercase() {
        Modification::try_from(raw_modification[0])
            .map_err(|_| format!("target_mod_code {value:?} has an invalid short code"))?
    } else if raw_modification.iter().all(u8::is_ascii_digit) {
        let raw = std::str::from_utf8(raw_modification)
            .map_err(|_| format!("target_mod_code {value:?} is not ASCII"))?;
        let id = raw
            .parse::<u32>()
            .map_err(|_| format!("target_mod_code {value:?} has an invalid ChEBI ID"))?;
        Modification::ChebiId(id)
    } else {
        return Err(format!(
            "target_mod_code {value:?} must select exactly one lowercase code or one numeric ChEBI ID"
        ));
    };

    let normalized = format_mod_code(unmodified_base, strand, modification);
    Ok(TargetModCode {
        unmodified_base,
        strand,
        modification,
        normalized,
    })
}

fn format_mod_code(
    unmodified_base: UnmodifiedBase,
    strand: MmStrand,
    modification: Modification,
) -> String {
    let base = char::from(u8::from(unmodified_base));
    let strand = match strand {
        MmStrand::Forward => '+',
        MmStrand::Reverse => '-',
    };
    let modification = match modification {
        Modification::Code(code) => char::from(code).to_string(),
        Modification::ChebiId(id) => id.to_string(),
    };
    format!("{base}{strand}{modification}")
}

fn validate_options(options: &ModBamOptions) -> Result<TargetModel, ModBamError> {
    crate::modification::types::validate_identifier("assay_id", &options.assay_id)
        .map_err(ModBamError::InvalidOptions)?;
    crate::modification::types::validate_identifier("sample", &options.sample)
        .map_err(ModBamError::InvalidOptions)?;
    if options.sample.contains(crate::sample::SAMPLE_DELIM) {
        return Err(ModBamError::InvalidOptions(format!(
            "sample {:?} must not contain {:?}",
            options.sample,
            crate::sample::SAMPLE_DELIM
        )));
    }
    let target =
        parse_target_mod_code(&options.target_mod_code).map_err(ModBamError::InvalidOptions)?;
    let candidate_rule = CandidateRule::parse(&options.candidate_rule, &target)
        .map_err(ModBamError::InvalidOptions)?;
    if let Some(threshold) = options.source_emission_threshold {
        if !threshold.is_finite() || !(0.0..=1.0).contains(&threshold) {
            return Err(ModBamError::InvalidOptions(format!(
                "source_emission_threshold must be finite and in [0, 1], got {threshold}"
            )));
        }
    }
    if options.mm_question_mark_policy == MmQuestionMarkPolicy::BelowEmissionThreshold
        && !matches!(options.source_emission_threshold, Some(threshold) if threshold > 0.0)
    {
        return Err(ModBamError::InvalidOptions(
            "mm_question_mark_policy=below_emission_threshold requires a positive source_emission_threshold"
                .to_owned(),
        ));
    }
    Ok(TargetModel {
        target,
        candidate_rule,
    })
}

#[derive(Clone, Debug)]
struct RecordIssue {
    reason: InvalidRecordReason,
    detail: String,
}

impl RecordIssue {
    fn new(reason: InvalidRecordReason, detail: impl Into<String>) -> Self {
        Self {
            reason,
            detail: detail.into(),
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct TargetCall {
    query_index: usize,
    ml_byte: u8,
    strand: MmStrand,
}

#[derive(Clone, Copy, Debug)]
struct TargetImplicitCall {
    query_index: usize,
    strand: MmStrand,
    state: ObservationState,
}

#[derive(Clone, Debug, Default)]
struct ParsedModTags {
    calls: Vec<TargetCall>,
    implicit_calls: Vec<TargetImplicitCall>,
    target_canonical_bases: usize,
    target_candidate_bases: usize,
    target_groups_skip_flag_omitted: usize,
    target_groups_low_probability: usize,
    target_groups_unknown: usize,
    implicit_low_probability_candidates: usize,
    unknown_candidates: usize,
    ml_values_consumed: usize,
}

fn validate_mm_delta_bounds(
    mm: &[u8],
    decoded_sequence: &[u8],
    is_reverse_complemented: bool,
) -> Result<(), RecordIssue> {
    for group in mm
        .split(|byte| *byte == b';')
        .filter(|group| !group.is_empty())
    {
        let Some(&canonical) = group.first() else {
            continue;
        };
        let stored_canonical = if is_reverse_complemented {
            match canonical {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' | b'U' => b'A',
                b'N' => b'N',
                _ => continue,
            }
        } else {
            canonical
        };
        let candidate_count = if canonical == b'N' {
            decoded_sequence.len()
        } else {
            decoded_sequence
                .iter()
                .filter(|&&base| base == stored_canonical)
                .count()
        };
        let Some(comma) = group.iter().position(|byte| *byte == b',') else {
            continue;
        };
        let mut canonical_index = 0usize;
        for raw_delta in group[comma + 1..].split(|byte| *byte == b',') {
            let delta_text = std::str::from_utf8(raw_delta).map_err(|_| {
                RecordIssue::new(InvalidRecordReason::InvalidMm, "MM delta is not ASCII")
            })?;
            let delta = delta_text.parse::<usize>().map_err(|_| {
                RecordIssue::new(
                    InvalidRecordReason::InvalidMm,
                    format!("invalid MM delta {delta_text:?}"),
                )
            })?;
            canonical_index = canonical_index.checked_add(delta).ok_or_else(|| {
                RecordIssue::new(
                    InvalidRecordReason::InvalidMm,
                    "MM delta index overflows usize",
                )
            })?;
            if canonical_index >= candidate_count {
                return Err(RecordIssue::new(
                    InvalidRecordReason::InvalidMm,
                    format!(
                        "MM delta refers beyond SEQ: canonical base {} call index {} but SEQ contains {candidate_count} candidates",
                        char::from(canonical),
                        canonical_index
                    ),
                ));
            }
            canonical_index += 1;
        }
    }
    Ok(())
}

fn parse_base_modifications(
    mm: &[u8],
    decoded_sequence: &[u8],
    is_reverse_complemented: bool,
) -> Result<BaseModifications, RecordIssue> {
    // noodles 0.85 recognizes `N` as an MM canonical base but decodes its
    // deltas by looking for literal `N` symbols in SEQ. SAM instead defines an
    // MM canonical `N` as counting every base. Parse one group at a time and
    // give only `N` groups an all-N coordinate model so noodles retains its
    // grammar/status/modification validation while producing the specified
    // any-base positions.
    let sequence = sam::alignment::record_buf::Sequence::from(decoded_sequence.to_vec());
    let any_base_sequence =
        sam::alignment::record_buf::Sequence::from(vec![b'N'; decoded_sequence.len()]);
    let mut groups: Vec<Group> = Vec::new();

    for raw_group in mm.split_inclusive(|byte| *byte == b';') {
        let group_sequence = if raw_group.first() == Some(&b'N') {
            &any_base_sequence
        } else {
            &sequence
        };
        let parsed =
            match BaseModifications::parse(raw_group, is_reverse_complemented, group_sequence) {
                Ok(parsed) => parsed,
                Err(error) => {
                    return Err(RecordIssue::new(
                        InvalidRecordReason::InvalidMm,
                        format!("cannot parse MM group: {error:?}"),
                    ));
                }
            };
        let mut parsed_groups: Vec<Group> = parsed.into();
        if parsed_groups.len() != 1 {
            return Err(RecordIssue::new(
                InvalidRecordReason::InvalidMm,
                "MM group parser did not return exactly one group",
            ));
        }
        groups.append(&mut parsed_groups);
    }

    Ok(BaseModifications::from(groups))
}

#[cfg(test)]
fn parse_mod_tags(
    mm: &[u8],
    ml: &[u8],
    mn: usize,
    decoded_sequence: &[u8],
    is_reverse_complemented: bool,
    mm_question_mark_policy: MmQuestionMarkPolicy,
    target: &TargetModCode,
) -> Result<ParsedModTags, RecordIssue> {
    let model = TargetModel {
        target: target.clone(),
        candidate_rule: CandidateRule::AllTargetCanonicalBases,
    };
    parse_mod_tags_with_model(
        mm,
        ml,
        mn,
        decoded_sequence,
        is_reverse_complemented,
        mm_question_mark_policy,
        &model,
    )
}

fn parse_mod_tags_with_model(
    mm: &[u8],
    ml: &[u8],
    mn: usize,
    decoded_sequence: &[u8],
    is_reverse_complemented: bool,
    mm_question_mark_policy: MmQuestionMarkPolicy,
    model: &TargetModel,
) -> Result<ParsedModTags, RecordIssue> {
    let target = &model.target;
    let candidate_rule = &model.candidate_rule;
    if mn != decoded_sequence.len() {
        return Err(RecordIssue::new(
            InvalidRecordReason::MnMismatch,
            format!(
                "MN is {mn}, but decoded SEQ length is {}",
                decoded_sequence.len()
            ),
        ));
    }
    validate_mm_delta_bounds(mm, decoded_sequence, is_reverse_complemented)?;

    // Construct from decoded bases. bam::record::Sequence::as_ref() exposes
    // packed nybbles and must never be passed to the MM parser.
    let base_modifications =
        parse_base_modifications(mm, decoded_sequence, is_reverse_complemented)?;

    let expected_ml_len = base_modifications
        .as_ref()
        .iter()
        .try_fold(0usize, |total, group| {
            let group_len = group
                .positions()
                .len()
                .checked_mul(group.modifications().len())
                .ok_or_else(|| {
                    RecordIssue::new(
                        InvalidRecordReason::MlLengthMismatch,
                        "MM group probability count overflows usize",
                    )
                })?;
            total.checked_add(group_len).ok_or_else(|| {
                RecordIssue::new(
                    InvalidRecordReason::MlLengthMismatch,
                    "total MM probability count overflows usize",
                )
            })
        })?;
    if ml.len() != expected_ml_len {
        return Err(RecordIssue::new(
            InvalidRecordReason::MlLengthMismatch,
            format!(
                "ML has {} values, but all MM groups require {expected_ml_len}",
                ml.len()
            ),
        ));
    }

    let target_canonical_bases = if target.unmodified_base == UnmodifiedBase::N {
        decoded_sequence.len()
    } else {
        let target_sequence_base = if is_reverse_complemented {
            u8::from(target.unmodified_base.complement())
        } else {
            u8::from(target.unmodified_base)
        };
        decoded_sequence
            .iter()
            .filter(|&&base| base == target_sequence_base)
            .count()
    };
    let target_candidate_positions =
        candidate_rule.candidate_positions(decoded_sequence, is_reverse_complemented, target);
    let target_candidate_set = target_candidate_positions
        .iter()
        .copied()
        .collect::<HashSet<_>>();
    let mut parsed = ParsedModTags {
        target_canonical_bases,
        target_candidate_bases: target_candidate_positions.len(),
        ml_values_consumed: expected_ml_len,
        ..ParsedModTags::default()
    };
    let mut ml_index = 0usize;
    let mut target_group_count = 0usize;

    for group in base_modifications.as_ref() {
        let group_sequence_base = if is_reverse_complemented {
            u8::from(group.unmodified_base().complement())
        } else {
            u8::from(group.unmodified_base())
        };
        for &query_index in group.positions() {
            let actual_base = decoded_sequence.get(query_index).ok_or_else(|| {
                RecordIssue::new(
                    InvalidRecordReason::InvalidMm,
                    format!(
                        "MM query index {query_index} exceeds decoded SEQ length {}",
                        decoded_sequence.len()
                    ),
                )
            })?;
            if group.unmodified_base() != UnmodifiedBase::N && *actual_base != group_sequence_base {
                return Err(RecordIssue::new(
                    InvalidRecordReason::InvalidMm,
                    format!(
                        "MM query index {query_index} points to base {}, expected canonical base {}",
                        char::from(*actual_base),
                        char::from(group_sequence_base)
                    ),
                ));
            }
        }
        let group_contains_target = group.unmodified_base() == target.unmodified_base
            && group.strand() == target.strand
            && group.modifications().contains(&target.modification);
        if group_contains_target {
            target_group_count += 1;
            if target_group_count > 1 {
                return Err(RecordIssue::new(
                    InvalidRecordReason::DuplicateTargetGroup,
                    format!(
                        "target code {} occurs in more than one MM group",
                        target.normalized
                    ),
                ));
            }
            if let Some(query_index) = group
                .positions()
                .iter()
                .find(|query_index| !target_candidate_set.contains(query_index))
            {
                return Err(RecordIssue::new(
                    InvalidRecordReason::CandidateRuleMismatch,
                    format!(
                        "explicit target call at query index {query_index} is outside candidate_rule {:?}",
                        candidate_rule.normalized()
                    ),
                ));
            }
            let omitted_count = target_candidate_positions
                .len()
                .checked_sub(group.positions().len())
                .ok_or_else(|| {
                    RecordIssue::new(
                        InvalidRecordReason::InvalidMm,
                        "MM lists more target positions than canonical bases",
                    )
                })?;
            let implicit_state = match group.status() {
                None => {
                    parsed.target_groups_skip_flag_omitted += 1;
                    parsed.implicit_low_probability_candidates += omitted_count;
                    ObservationState::ImplicitBelowEmissionThreshold
                }
                Some(Status::Implicit) => {
                    parsed.target_groups_low_probability += 1;
                    parsed.implicit_low_probability_candidates += omitted_count;
                    ObservationState::ImplicitBelowEmissionThreshold
                }
                Some(Status::Explicit) => {
                    parsed.target_groups_unknown += 1;
                    match mm_question_mark_policy {
                        MmQuestionMarkPolicy::Unknown => {
                            parsed.unknown_candidates += omitted_count;
                            ObservationState::Unknown
                        }
                        MmQuestionMarkPolicy::BelowEmissionThreshold => {
                            parsed.implicit_low_probability_candidates += omitted_count;
                            ObservationState::ImplicitBelowEmissionThreshold
                        }
                    }
                }
            };
            let listed_positions = group.positions().iter().copied().collect::<HashSet<_>>();
            parsed.implicit_calls.extend(
                target_candidate_positions
                    .iter()
                    .copied()
                    .filter(|query_index| !listed_positions.contains(query_index))
                    .map(|query_index| TargetImplicitCall {
                        query_index,
                        strand: group.strand(),
                        state: implicit_state,
                    }),
            );
            if parsed.implicit_calls.len()
                != parsed.implicit_low_probability_candidates + parsed.unknown_candidates
            {
                return Err(RecordIssue::new(
                    InvalidRecordReason::InvalidMm,
                    "MM omitted-candidate accounting is inconsistent",
                ));
            }
        }

        // SAM ML ordering is group-major, then position-major, then code-major.
        // Advance for every code before deciding whether this is the target.
        for &query_index in group.positions() {
            for &modification in group.modifications() {
                let ml_byte = ml[ml_index];
                ml_index += 1;
                if group.unmodified_base() == target.unmodified_base
                    && group.strand() == target.strand
                    && modification == target.modification
                {
                    parsed.calls.push(TargetCall {
                        query_index,
                        ml_byte,
                        strand: group.strand(),
                    });
                }
            }
        }
    }
    debug_assert_eq!(ml_index, ml.len());

    let mut seen_target_positions = HashSet::with_capacity(parsed.calls.len());
    if let Some(duplicate) = parsed
        .calls
        .iter()
        .map(|call| call.query_index)
        .find(|&query_index| !seen_target_positions.insert(query_index))
    {
        return Err(RecordIssue::new(
            InvalidRecordReason::DuplicateTargetObservation,
            format!("target code occurs more than once at query index {duplicate}"),
        ));
    }

    Ok(parsed)
}

#[derive(Debug, Default)]
struct RawModTags {
    mm: Option<Vec<u8>>,
    ml: Option<Vec<u8>>,
    mn: Option<usize>,
}

fn decode_raw_mod_tags(record: &bam::Record) -> Result<RawModTags, RecordIssue> {
    let mut tags = RawModTags::default();

    for result in record.data().iter() {
        let (tag, value) = result.map_err(|error| {
            RecordIssue::new(
                InvalidRecordReason::InvalidAuxData,
                format!("cannot decode auxiliary field: {error}"),
            )
        })?;

        match tag {
            Tag::BASE_MODIFICATIONS => {
                if tags.mm.is_some() {
                    return Err(RecordIssue::new(
                        InvalidRecordReason::DuplicateMm,
                        "MM occurs more than once",
                    ));
                }
                let Value::String(value) = value else {
                    return Err(RecordIssue::new(
                        InvalidRecordReason::InvalidMmType,
                        format!("MM must be Z, found {:?}", value.ty()),
                    ));
                };
                tags.mm = Some(value.iter().copied().collect());
            }
            Tag::BASE_MODIFICATION_PROBABILITIES => {
                if tags.ml.is_some() {
                    return Err(RecordIssue::new(
                        InvalidRecordReason::DuplicateMl,
                        "ML occurs more than once",
                    ));
                }
                let Value::Array(Array::UInt8(values)) = value else {
                    return Err(RecordIssue::new(
                        InvalidRecordReason::InvalidMlType,
                        "ML must be B:C",
                    ));
                };
                let values = values
                    .iter()
                    .collect::<io::Result<Vec<_>>>()
                    .map_err(|error| {
                        RecordIssue::new(
                            InvalidRecordReason::InvalidMlValue,
                            format!("cannot decode ML value: {error}"),
                        )
                    })?;
                tags.ml = Some(values);
            }
            Tag::BASE_MODIFICATION_SEQUENCE_LENGTH => {
                if tags.mn.is_some() {
                    return Err(RecordIssue::new(
                        InvalidRecordReason::DuplicateMn,
                        "MN occurs more than once",
                    ));
                }
                let value = value.as_int().ok_or_else(|| {
                    RecordIssue::new(InvalidRecordReason::InvalidMn, "MN must be an integer")
                })?;
                let value = usize::try_from(value).map_err(|_| {
                    RecordIssue::new(
                        InvalidRecordReason::InvalidMn,
                        format!("MN must be nonnegative and fit usize, got {value}"),
                    )
                })?;
                tags.mn = Some(value);
            }
            _ => {}
        }
    }

    if tags.mm.is_none() {
        return Err(RecordIssue::new(
            InvalidRecordReason::MissingMm,
            "required MM tag is absent",
        ));
    }
    if tags.ml.is_none() {
        return Err(RecordIssue::new(
            InvalidRecordReason::MissingMl,
            "required ML tag is absent",
        ));
    }
    if tags.mn.is_none() {
        return Err(RecordIssue::new(
            InvalidRecordReason::MissingMn,
            "required MN tag is absent",
        ));
    }
    Ok(tags)
}

#[derive(Debug, Default)]
struct RecordImport {
    observations: Vec<ModObservation>,
    target_group_seen: bool,
    target_canonical_bases: usize,
    target_candidate_bases: usize,
    target_groups_skip_flag_omitted: usize,
    target_groups_low_probability: usize,
    target_groups_unknown: usize,
    implicit_low_probability_candidates: usize,
    unknown_candidates: usize,
    ml_values_consumed: usize,
    explicit_target_calls: usize,
    target_calls_in_insertions: usize,
    target_calls_in_soft_clips: usize,
    implicit_calls_in_insertions: usize,
    implicit_calls_in_soft_clips: usize,
}

fn genomic_strand(is_reverse_complemented: bool, mm_strand: MmStrand) -> Strand {
    match (is_reverse_complemented, mm_strand) {
        (false, MmStrand::Forward) | (true, MmStrand::Reverse) => Strand::Plus,
        (false, MmStrand::Reverse) | (true, MmStrand::Forward) => Strand::Minus,
    }
}

#[cfg(test)]
fn normalize_target_calls(
    chrom: &str,
    read_name: &str,
    is_reverse_complemented: bool,
    projection: &QueryReferenceProjection,
    parsed: ParsedModTags,
    options: &ModBamOptions,
    target: &TargetModCode,
) -> Result<RecordImport, RecordIssue> {
    let model = TargetModel {
        target: target.clone(),
        candidate_rule: CandidateRule::AllTargetCanonicalBases,
    };
    normalize_target_calls_with_model(
        chrom,
        read_name,
        is_reverse_complemented,
        projection,
        parsed,
        options,
        &model,
    )
}

fn normalize_target_calls_with_model(
    chrom: &str,
    read_name: &str,
    is_reverse_complemented: bool,
    projection: &QueryReferenceProjection,
    parsed: ParsedModTags,
    options: &ModBamOptions,
    model: &TargetModel,
) -> Result<RecordImport, RecordIssue> {
    let target = &model.target;
    let candidate_rule = &model.candidate_rule;
    let target_group_seen = parsed.target_groups_skip_flag_omitted
        + parsed.target_groups_low_probability
        + parsed.target_groups_unknown
        > 0;
    let mut imported = RecordImport {
        target_group_seen,
        target_canonical_bases: parsed.target_canonical_bases,
        target_candidate_bases: parsed.target_candidate_bases,
        target_groups_skip_flag_omitted: parsed.target_groups_skip_flag_omitted,
        target_groups_low_probability: parsed.target_groups_low_probability,
        target_groups_unknown: parsed.target_groups_unknown,
        implicit_low_probability_candidates: parsed.implicit_low_probability_candidates,
        unknown_candidates: parsed.unknown_candidates,
        ml_values_consumed: parsed.ml_values_consumed,
        explicit_target_calls: parsed.calls.len(),
        ..RecordImport::default()
    };
    let read_id = if read_name.contains(crate::sample::SAMPLE_DELIM) {
        let (sample, _) = crate::sample::split_tagged_read_name(read_name).ok_or_else(|| {
            RecordIssue::new(
                InvalidRecordReason::InvalidReadName,
                format!(
                    "query name {read_name:?} has a malformed sample prefix using {:?}",
                    crate::sample::SAMPLE_DELIM
                ),
            )
        })?;
        if sample != options.sample {
            return Err(RecordIssue::new(
                InvalidRecordReason::InvalidReadName,
                format!(
                    "query name {read_name:?} has sample prefix {sample:?}, expected {:?}",
                    options.sample
                ),
            ));
        }
        read_name.to_owned()
    } else {
        crate::sample::tagged_read_name(&options.sample, read_name)
    };

    for call in parsed.calls {
        let projected = projection
            .query_bases
            .get(call.query_index)
            .ok_or_else(|| {
                RecordIssue::new(
                    InvalidRecordReason::InvalidProjection,
                    format!(
                        "MM query index {} exceeds projected SEQ length {}",
                        call.query_index,
                        projection.query_bases.len()
                    ),
                )
            })?;
        let pos0 = match projected {
            QueryBaseProjection::Reference(pos0) => *pos0,
            QueryBaseProjection::Insertion => {
                imported.target_calls_in_insertions += 1;
                continue;
            }
            QueryBaseProjection::SoftClip => {
                imported.target_calls_in_soft_clips += 1;
                continue;
            }
        };

        let observation = ModObservation {
            key: ModObservationKey {
                assay_id: options.assay_id.clone(),
                sample: options.sample.clone(),
                read_id: read_id.clone(),
                site: ModSiteKey {
                    chrom: chrom.to_owned(),
                    pos0,
                    strand: genomic_strand(is_reverse_complemented, call.strand),
                    mod_code: target.normalized.clone(),
                },
            },
            probability: Some(ml_byte_to_probability(call.ml_byte)),
            observation_state: ObservationState::ExplicitProbability,
            context: candidate_rule.observation_context().map(str::to_owned),
            source_transcript_id: None,
            source_pos0: None,
        };
        observation
            .validate()
            .map_err(|error| RecordIssue::new(InvalidRecordReason::InvalidObservation, error))?;
        imported.observations.push(observation);
    }

    for call in parsed.implicit_calls {
        let projected = projection
            .query_bases
            .get(call.query_index)
            .ok_or_else(|| {
                RecordIssue::new(
                    InvalidRecordReason::InvalidProjection,
                    format!(
                        "MM implicit query index {} exceeds projected SEQ length {}",
                        call.query_index,
                        projection.query_bases.len()
                    ),
                )
            })?;
        let pos0 = match projected {
            QueryBaseProjection::Reference(pos0) => *pos0,
            QueryBaseProjection::Insertion => {
                imported.implicit_calls_in_insertions += 1;
                continue;
            }
            QueryBaseProjection::SoftClip => {
                imported.implicit_calls_in_soft_clips += 1;
                continue;
            }
        };
        let observation = ModObservation {
            key: ModObservationKey {
                assay_id: options.assay_id.clone(),
                sample: options.sample.clone(),
                read_id: read_id.clone(),
                site: ModSiteKey {
                    chrom: chrom.to_owned(),
                    pos0,
                    strand: genomic_strand(is_reverse_complemented, call.strand),
                    mod_code: target.normalized.clone(),
                },
            },
            probability: None,
            observation_state: call.state,
            context: candidate_rule.observation_context().map(str::to_owned),
            source_transcript_id: None,
            source_pos0: None,
        };
        observation
            .validate()
            .map_err(|error| RecordIssue::new(InvalidRecordReason::InvalidObservation, error))?;
        imported.observations.push(observation);
    }

    imported
        .observations
        .sort_by(|left, right| left.key.cmp(&right.key));
    if let Some(pair) = imported
        .observations
        .windows(2)
        .find(|pair| pair[0].key == pair[1].key)
    {
        return Err(RecordIssue::new(
            InvalidRecordReason::DuplicateTargetObservation,
            format!("duplicate target observation key {:?}", pair[0].key),
        ));
    }

    Ok(imported)
}

fn import_record(
    record: &bam::Record,
    header: &sam::Header,
    read_name: &str,
    options: &ModBamOptions,
    model: &TargetModel,
) -> Result<RecordImport, RecordIssue> {
    let reference_id = record
        .reference_sequence_id()
        .transpose()
        .map_err(|error| {
            RecordIssue::new(
                InvalidRecordReason::InvalidReference,
                format!("cannot decode reference ID: {error}"),
            )
        })?
        .ok_or_else(|| {
            RecordIssue::new(
                InvalidRecordReason::InvalidReference,
                "mapped record has no reference ID",
            )
        })?;
    let (reference_name, reference_sequence) = header
        .reference_sequences()
        .get_index(reference_id)
        .ok_or_else(|| {
            RecordIssue::new(
                InvalidRecordReason::InvalidReference,
                format!("reference ID {reference_id} is absent from BAM header"),
            )
        })?;
    let chrom = String::from_utf8(reference_name.iter().copied().collect()).map_err(|_| {
        RecordIssue::new(
            InvalidRecordReason::InvalidReference,
            "reference name is not valid UTF-8",
        )
    })?;
    if chrom.trim().is_empty() || chrom == "*" || chrom.chars().any(char::is_control) {
        return Err(RecordIssue::new(
            InvalidRecordReason::InvalidReference,
            format!("invalid reference name {chrom:?}"),
        ));
    }
    let reference_len = usize::from(reference_sequence.length());

    let alignment_start = record
        .alignment_start()
        .transpose()
        .map_err(|error| {
            RecordIssue::new(
                InvalidRecordReason::InvalidAlignmentStart,
                format!("cannot decode alignment start: {error}"),
            )
        })?
        .ok_or_else(|| {
            RecordIssue::new(
                InvalidRecordReason::InvalidAlignmentStart,
                "mapped record has no alignment start",
            )
        })?;
    let alignment_start0 = u32::try_from(usize::from(alignment_start) - 1).map_err(|_| {
        RecordIssue::new(
            InvalidRecordReason::InvalidAlignmentStart,
            "alignment start exceeds u32",
        )
    })?;

    let decoded_sequence = record.sequence().iter().collect::<Vec<_>>();
    let raw_tags = decode_raw_mod_tags(record)?;
    let parsed = parse_mod_tags_with_model(
        raw_tags.mm.as_deref().expect("validated MM"),
        raw_tags.ml.as_deref().expect("validated ML"),
        raw_tags.mn.expect("validated MN"),
        &decoded_sequence,
        record.flags().is_reverse_complemented(),
        options.mm_question_mark_policy,
        model,
    )?;

    let cigar = record
        .cigar()
        .iter()
        .collect::<io::Result<Vec<_>>>()
        .map_err(|error| {
            RecordIssue::new(
                InvalidRecordReason::InvalidCigar,
                format!("cannot decode CIGAR: {error}"),
            )
        })?;
    let projection = project_query_to_reference(alignment_start0, decoded_sequence.len(), &cigar)
        .map_err(|error| {
        RecordIssue::new(InvalidRecordReason::InvalidProjection, error.to_string())
    })?;
    if usize::try_from(projection.reference_end0).unwrap_or(usize::MAX) > reference_len {
        return Err(RecordIssue::new(
            InvalidRecordReason::ReferenceOutOfBounds,
            format!(
                "CIGAR ends at {} beyond reference {chrom:?} length {reference_len}",
                projection.reference_end0
            ),
        ));
    }

    normalize_target_calls_with_model(
        &chrom,
        read_name,
        record.flags().is_reverse_complemented(),
        &projection,
        parsed,
        options,
        model,
    )
}

#[derive(Debug)]
enum PrimaryOutcome {
    BelowMapq,
    Invalid(RecordIssue),
    Valid(RecordImport),
}

#[derive(Debug)]
struct PrimaryEntry {
    first_ordinal: usize,
    count: usize,
    outcome: PrimaryOutcome,
}

fn fail_record(record_ordinal: usize, read_name: Option<&str>, issue: RecordIssue) -> ModBamError {
    ModBamError::InvalidRecord {
        record_ordinal,
        read_name: read_name.unwrap_or("<unknown>").to_owned(),
        reason: issue.reason,
        detail: issue.detail,
    }
}

fn merge_record_import(
    qc: &mut ModBamQc,
    observations: &mut Vec<ModObservation>,
    record: RecordImport,
) {
    qc.retained_records += 1;
    qc.target_canonical_bases += record.target_canonical_bases;
    qc.target_candidate_bases += record.target_candidate_bases;
    if !record.target_group_seen && record.target_candidate_bases > 0 {
        qc.records_without_target_group += 1;
    }
    qc.target_groups_skip_flag_omitted += record.target_groups_skip_flag_omitted;
    qc.target_groups_low_probability += record.target_groups_low_probability;
    qc.target_groups_unknown += record.target_groups_unknown;
    qc.implicit_low_probability_candidates += record.implicit_low_probability_candidates;
    qc.unknown_candidates += record.unknown_candidates;
    qc.ml_values_consumed += record.ml_values_consumed;
    qc.explicit_target_calls += record.explicit_target_calls;
    qc.target_calls_in_insertions += record.target_calls_in_insertions;
    qc.target_calls_in_soft_clips += record.target_calls_in_soft_clips;
    qc.implicit_calls_in_insertions += record.implicit_calls_in_insertions;
    qc.implicit_calls_in_soft_clips += record.implicit_calls_in_soft_clips;
    observations.extend(record.observations);
}

fn scan_modbam<R: io::Read>(
    source: R,
    path: &Path,
    options: &ModBamOptions,
    model: &TargetModel,
) -> Result<ModBamImportResult, ModBamError> {
    let mut reader = bam::io::Reader::new(source);
    let header = reader
        .read_header()
        .map_err(|source| ModBamError::ReadHeader {
            path: path.to_owned(),
            source,
        })?;
    let mut qc = ModBamQc::default();
    let mut primary_by_name: HashMap<String, PrimaryEntry> = HashMap::new();

    for (record_index, result) in reader.records().enumerate() {
        let record_ordinal = record_index + 1;
        qc.total_records += 1;
        let record = match result {
            Ok(record) => record,
            Err(error) => {
                let issue = RecordIssue::new(
                    InvalidRecordReason::RecordDecode,
                    format!("cannot decode BAM record: {error}"),
                );
                return Err(fail_record(record_ordinal, None, issue));
            }
        };

        let flags = record.flags();
        if flags.is_unmapped() {
            qc.skipped_unmapped += 1;
            continue;
        }
        if flags.is_secondary() {
            qc.skipped_secondary += 1;
            continue;
        }
        if flags.is_supplementary() {
            qc.skipped_supplementary += 1;
            continue;
        }

        let read_name = match record
            .name()
            .map(|name| String::from_utf8(name.iter().copied().collect::<Vec<_>>()))
        {
            Some(Err(_)) => {
                let issue = RecordIssue::new(
                    InvalidRecordReason::InvalidReadName,
                    "query name is not valid UTF-8",
                );
                if options.invalid_record_policy == InvalidRecordPolicy::Fail {
                    return Err(fail_record(record_ordinal, None, issue));
                }
                *qc.invalid_records.entry(issue.reason).or_default() += 1;
                continue;
            }
            Some(Ok(name)) if name.trim().is_empty() || name == "*" => {
                let issue = RecordIssue::new(
                    InvalidRecordReason::MissingReadName,
                    "mapped primary record has no query name",
                );
                if options.invalid_record_policy == InvalidRecordPolicy::Fail {
                    return Err(fail_record(record_ordinal, None, issue));
                }
                *qc.invalid_records.entry(issue.reason).or_default() += 1;
                continue;
            }
            Some(Ok(name)) => {
                if name.chars().any(char::is_control) {
                    let issue = RecordIssue::new(
                        InvalidRecordReason::InvalidReadName,
                        format!("query name contains a control character: {name:?}"),
                    );
                    if options.invalid_record_policy == InvalidRecordPolicy::Fail {
                        return Err(fail_record(record_ordinal, Some(&name), issue));
                    }
                    *qc.invalid_records.entry(issue.reason).or_default() += 1;
                    continue;
                } else {
                    name
                }
            }
            None => {
                let issue = RecordIssue::new(
                    InvalidRecordReason::MissingReadName,
                    "mapped primary record has no query name",
                );
                if options.invalid_record_policy == InvalidRecordPolicy::Fail {
                    return Err(fail_record(record_ordinal, None, issue));
                }
                *qc.invalid_records.entry(issue.reason).or_default() += 1;
                continue;
            }
        };

        if let Some(entry) = primary_by_name.get_mut(&read_name) {
            entry.count += 1;
            if options.invalid_record_policy == InvalidRecordPolicy::Fail {
                return Err(fail_record(
                    record_ordinal,
                    Some(&read_name),
                    RecordIssue::new(
                        InvalidRecordReason::DuplicatePrimary,
                        format!(
                            "primary alignment first occurred at record {}",
                            entry.first_ordinal
                        ),
                    ),
                ));
            }
            continue;
        }

        let mapq = record.mapping_quality().map(u8::from).unwrap_or(0);
        let outcome = if mapq < options.min_mapq {
            PrimaryOutcome::BelowMapq
        } else {
            match import_record(&record, &header, &read_name, options, model) {
                Ok(imported) => PrimaryOutcome::Valid(imported),
                Err(issue) if options.invalid_record_policy == InvalidRecordPolicy::Skip => {
                    PrimaryOutcome::Invalid(issue)
                }
                Err(issue) => return Err(fail_record(record_ordinal, Some(&read_name), issue)),
            }
        };
        primary_by_name.insert(
            read_name,
            PrimaryEntry {
                first_ordinal: record_ordinal,
                count: 1,
                outcome,
            },
        );
    }

    let mut observations = Vec::new();
    for entry in primary_by_name.into_values() {
        if entry.count > 1 {
            qc.duplicate_primary_reads += 1;
            qc.skipped_duplicate_primary_records += entry.count;
            *qc.invalid_records
                .entry(InvalidRecordReason::DuplicatePrimary)
                .or_default() += entry.count;
            continue;
        }
        match entry.outcome {
            PrimaryOutcome::BelowMapq => qc.skipped_below_mapq += 1,
            PrimaryOutcome::Invalid(issue) => {
                *qc.invalid_records.entry(issue.reason).or_default() += 1;
            }
            PrimaryOutcome::Valid(record) => {
                merge_record_import(&mut qc, &mut observations, record)
            }
        }
    }
    observations.sort_by(|left, right| left.key.cmp(&right.key));
    qc.emitted_observations = observations.len();
    for observation in &observations {
        match observation.observation_state {
            ObservationState::ExplicitProbability => qc.emitted_explicit_observations += 1,
            ObservationState::ImplicitBelowEmissionThreshold => {
                qc.emitted_implicit_observations += 1
            }
            ObservationState::Unknown => qc.emitted_unknown_observations += 1,
        }
    }

    let omitted_candidate_count = qc.implicit_low_probability_candidates + qc.unknown_candidates;
    let implicit_skip_policy = if qc.unknown_candidates > 0 {
        ImplicitSkipPolicy::Unknown
    } else if qc.implicit_low_probability_candidates > 0
        || options.mm_question_mark_policy == MmQuestionMarkPolicy::BelowEmissionThreshold
    {
        ImplicitSkipPolicy::LowProbability
    } else {
        ImplicitSkipPolicy::NotApplicable
    };
    let candidate_observations_complete = qc.retained_records > 0
        && qc.records_without_target_group == 0
        && qc.invalid_records.is_empty();

    Ok(ModBamImportResult {
        observations,
        qc,
        semantics: ModBamImportSemantics {
            target_mod_code: model.target.normalized.clone(),
            candidate_rule: model.candidate_rule.normalized().to_owned(),
            source_emission_threshold: options.source_emission_threshold,
            mm_question_mark_policy: options.mm_question_mark_policy,
            implicit_skip_policy,
            candidate_observations_complete,
            explicit_observations_only: omitted_candidate_count == 0,
            ml_probability_semantics: MlProbabilitySemantics::SamIntervalMidpoint,
        },
    })
}

/// Import a genome-aligned Dorado/modBAM.
///
/// Unmapped, secondary, and supplementary records are always excluded in V1.
/// Duplicate mapped primary query names are fatal under fail policy; under skip
/// policy, every primary record for that query name is discarded.
pub fn read_modbam(
    path: &Path,
    options: &ModBamOptions,
) -> Result<ModBamImportResult, ModBamError> {
    let model = validate_options(options)?;
    let file = File::open(path).map_err(|source| ModBamError::Open {
        path: path.to_owned(),
        source,
    })?;
    scan_modbam(file, path, options, &model)
}

#[cfg(test)]
mod tests {
    use std::num::NonZero;
    use std::time::{SystemTime, UNIX_EPOCH};

    use sam::alignment::io::Write as _;
    use sam::alignment::record::Flags;

    use super::*;

    fn target(value: &str) -> TargetModCode {
        parse_target_mod_code(value).unwrap()
    }

    fn options() -> ModBamOptions {
        ModBamOptions::new("rna004_m6a", "S1", "A+a")
    }

    #[test]
    fn ml_probability_uses_sam_interval_midpoint() {
        assert_eq!(ml_probability_interval(0), (0.0, 1.0 / 256.0));
        assert_eq!(ml_probability_interval(255), (255.0 / 256.0, 1.0));
        assert_eq!(ml_byte_to_probability(0), 0.5 / 256.0);
        assert_eq!(ml_byte_to_probability(255), 255.5 / 256.0);
        assert_eq!(
            MlProbabilitySemantics::SamIntervalMidpoint.as_str(),
            "sam_ml_interval_midpoint_v1"
        );
        assert_eq!(MmQuestionMarkPolicy::Unknown.as_str(), "unknown");
        assert_eq!(
            MmQuestionMarkPolicy::BelowEmissionThreshold.as_str(),
            "below_emission_threshold"
        );
    }

    #[test]
    fn question_mark_threshold_override_requires_a_positive_threshold() {
        let mut options = options();
        options.mm_question_mark_policy = MmQuestionMarkPolicy::BelowEmissionThreshold;
        assert!(matches!(
            validate_options(&options),
            Err(ModBamError::InvalidOptions(message))
                if message.contains("positive source_emission_threshold")
        ));

        options.source_emission_threshold = Some(0.0);
        assert!(validate_options(&options).is_err());

        options.source_emission_threshold = Some(0.05);
        assert!(validate_options(&options).is_ok());
    }

    #[test]
    fn centered_iupac_candidate_rules_are_validated_and_normalized() {
        let mut options = options();
        options.candidate_rule = "drach".to_owned();
        let model = validate_options(&options).unwrap();
        assert_eq!(model.candidate_rule.normalized(), "DRACH");

        for invalid in ["RA", "DRCCH", "DRA?H", "NNCNN"] {
            options.candidate_rule = invalid.to_owned();
            assert!(
                validate_options(&options).is_err(),
                "{invalid:?} unexpectedly passed"
            );
        }
    }

    #[test]
    fn motif_rule_expands_only_matching_candidates_in_both_bam_orientations() {
        let target = target("A+a");
        let model = TargetModel {
            candidate_rule: CandidateRule::parse("DRACH", &target).unwrap(),
            target,
        };
        let forward = parse_mod_tags_with_model(
            b"A+a.;",
            &[],
            5,
            b"AAACA",
            false,
            MmQuestionMarkPolicy::Unknown,
            &model,
        )
        .unwrap();
        assert_eq!(forward.target_canonical_bases, 4);
        assert_eq!(forward.target_candidate_bases, 1);
        assert_eq!(forward.implicit_calls.len(), 1);
        assert_eq!(forward.implicit_calls[0].query_index, 2);
        let projection = project_query_to_reference(100, 5, &[Op::new(Kind::Match, 5)]).unwrap();
        let imported = normalize_target_calls_with_model(
            "chr1",
            "motif",
            false,
            &projection,
            forward,
            &options(),
            &model,
        )
        .unwrap();
        assert_eq!(imported.observations.len(), 1);
        assert_eq!(imported.observations[0].context.as_deref(), Some("DRACH"));

        let reverse = parse_mod_tags_with_model(
            b"A+a.;",
            &[],
            5,
            b"TGTTT",
            true,
            MmQuestionMarkPolicy::Unknown,
            &model,
        )
        .unwrap();
        assert_eq!(reverse.target_canonical_bases, 4);
        assert_eq!(reverse.target_candidate_bases, 1);
        assert_eq!(reverse.implicit_calls[0].query_index, 2);
    }

    #[test]
    fn motif_rule_rejects_explicit_calls_outside_declared_candidates() {
        let target = target("A+a");
        let model = TargetModel {
            candidate_rule: CandidateRule::parse("DRACH", &target).unwrap(),
            target,
        };
        let error = parse_mod_tags_with_model(
            b"A+a.,0;",
            &[255],
            5,
            b"AAACA",
            false,
            MmQuestionMarkPolicy::Unknown,
            &model,
        )
        .unwrap_err();
        assert!(matches!(
            error.reason,
            InvalidRecordReason::CandidateRuleMismatch
        ));
    }

    #[test]
    fn cigar_projects_splice_insertion_soft_clip_and_all_kinds() {
        let cigar = [
            Op::new(Kind::SoftClip, 2),
            Op::new(Kind::Match, 3),
            Op::new(Kind::Insertion, 2),
            Op::new(Kind::SequenceMatch, 1),
            Op::new(Kind::Skip, 4),
            Op::new(Kind::SequenceMismatch, 2),
            Op::new(Kind::Deletion, 1),
            Op::new(Kind::Match, 1),
            Op::new(Kind::HardClip, 1),
            Op::new(Kind::Pad, 1),
        ];
        let actual = project_query_to_reference(100, 11, &cigar).unwrap();
        assert_eq!(
            actual.query_bases,
            vec![
                QueryBaseProjection::SoftClip,
                QueryBaseProjection::SoftClip,
                QueryBaseProjection::Reference(100),
                QueryBaseProjection::Reference(101),
                QueryBaseProjection::Reference(102),
                QueryBaseProjection::Insertion,
                QueryBaseProjection::Insertion,
                QueryBaseProjection::Reference(103),
                QueryBaseProjection::Reference(108),
                QueryBaseProjection::Reference(109),
                QueryBaseProjection::Reference(111),
            ]
        );
        assert_eq!(actual.reference_end0, 112);
    }

    #[test]
    fn cigar_rejects_query_length_mismatch_and_zero_length() {
        assert!(matches!(
            project_query_to_reference(0, 3, &[Op::new(Kind::Match, 2)]),
            Err(ProjectionError::QueryLengthMismatch { .. })
        ));
        assert!(matches!(
            project_query_to_reference(0, 1, &[Op::new(Kind::Match, 0)]),
            Err(ProjectionError::ZeroLength { .. })
        ));
    }

    #[test]
    fn parser_consumes_multi_mod_ml_before_target_filtering() {
        let parsed = parse_mod_tags(
            b"C+mh,0,1;A+a.,0;",
            &[10, 20, 30, 40, 50],
            5,
            b"CACCA",
            false,
            MmQuestionMarkPolicy::Unknown,
            &target("C+h"),
        )
        .unwrap();

        assert_eq!(parsed.ml_values_consumed, 5);
        assert_eq!(parsed.calls.len(), 2);
        assert_eq!(parsed.calls[0].query_index, 0);
        assert_eq!(parsed.calls[0].ml_byte, 20);
        assert_eq!(parsed.calls[1].query_index, 3);
        assert_eq!(parsed.calls[1].ml_byte, 40);
        assert_eq!(parsed.target_groups_skip_flag_omitted, 1);
    }

    #[test]
    fn parser_uses_decoded_sequence_and_flag_orientation() {
        let forward = parse_mod_tags(
            b"A+a,1;",
            &[200],
            5,
            b"AACAA",
            false,
            MmQuestionMarkPolicy::Unknown,
            &target("A+a"),
        )
        .unwrap();
        assert_eq!(forward.calls[0].query_index, 1);

        // Stored BAM SEQ is the reverse complement of the as-sequenced bases.
        let reverse = parse_mod_tags(
            b"A+a,1;",
            &[200],
            5,
            b"TTGTT",
            true,
            MmQuestionMarkPolicy::Unknown,
            &target("A+a"),
        )
        .unwrap();
        assert_eq!(reverse.calls[0].query_index, 3);
    }

    #[test]
    fn parser_treats_n_canonical_groups_as_any_base_in_both_orientations() {
        let forward = parse_mod_tags(
            b"N+n.,1;",
            &[200],
            4,
            b"ACGT",
            false,
            MmQuestionMarkPolicy::Unknown,
            &target("N+n"),
        )
        .unwrap();
        assert_eq!(forward.target_canonical_bases, 4);
        assert_eq!(forward.target_candidate_bases, 4);
        assert_eq!(forward.calls[0].query_index, 1);
        assert_eq!(forward.implicit_calls.len(), 3);

        // An N group counts from the original 5' end even when stored SEQ is
        // reverse complemented, so delta 1 resolves to stored query index 2.
        let reverse = parse_mod_tags(
            b"N+n.,1;",
            &[200],
            4,
            b"ACGT",
            true,
            MmQuestionMarkPolicy::Unknown,
            &target("N+n"),
        )
        .unwrap();
        assert_eq!(reverse.target_canonical_bases, 4);
        assert_eq!(reverse.target_candidate_bases, 4);
        assert_eq!(reverse.calls[0].query_index, 2);
        assert_eq!(reverse.implicit_calls.len(), 3);
    }

    #[test]
    fn non_target_n_group_does_not_block_target_import_or_ml_ordering() {
        let parsed = parse_mod_tags(
            b"A+a.,0;N+n.,0;",
            &[255, 200],
            4,
            b"ACGT",
            false,
            MmQuestionMarkPolicy::Unknown,
            &target("A+a"),
        )
        .unwrap();

        assert_eq!(parsed.ml_values_consumed, 2);
        assert_eq!(parsed.calls.len(), 1);
        assert_eq!(parsed.calls[0].query_index, 0);
        assert_eq!(parsed.calls[0].ml_byte, 255);
    }

    #[test]
    fn reverse_flag_call_projects_without_a_second_reversal() {
        let target = target("A+a");
        let parsed = parse_mod_tags(
            b"A+a,1;",
            &[200],
            5,
            b"TTGTT",
            true,
            MmQuestionMarkPolicy::Unknown,
            &target,
        )
        .unwrap();
        let projection = project_query_to_reference(100, 5, &[Op::new(Kind::Match, 5)]).unwrap();
        let imported = normalize_target_calls(
            "chr1",
            "reverse",
            true,
            &projection,
            parsed,
            &options(),
            &target,
        )
        .unwrap();

        assert_eq!(imported.observations.len(), 4);
        let explicit = imported
            .observations
            .iter()
            .find(|observation| {
                observation.observation_state == ObservationState::ExplicitProbability
            })
            .unwrap();
        assert_eq!(explicit.key.site.pos0, 103);
        assert_eq!(explicit.key.site.strand, Strand::Minus);
    }

    #[test]
    fn parser_rejects_mn_and_ml_mismatches() {
        let mn_error = parse_mod_tags(
            b"A+a,0;",
            &[1],
            3,
            b"AAAA",
            false,
            MmQuestionMarkPolicy::Unknown,
            &target("A+a"),
        )
        .unwrap_err();
        assert_eq!(mn_error.reason, InvalidRecordReason::MnMismatch);

        let ml_error = parse_mod_tags(
            b"A+ah,0,0;",
            &[1, 2, 3],
            4,
            b"AAAA",
            false,
            MmQuestionMarkPolicy::Unknown,
            &target("A+a"),
        )
        .unwrap_err();
        assert_eq!(ml_error.reason, InvalidRecordReason::MlLengthMismatch);

        let duplicate_group = parse_mod_tags(
            b"A+a,0;A+a,1;",
            &[1, 2],
            2,
            b"AA",
            false,
            MmQuestionMarkPolicy::Unknown,
            &target("A+a"),
        )
        .unwrap_err();
        assert_eq!(
            duplicate_group.reason,
            InvalidRecordReason::DuplicateTargetGroup
        );

        let beyond = parse_mod_tags(
            b"A+a,1;",
            &[1],
            1,
            b"A",
            false,
            MmQuestionMarkPolicy::Unknown,
            &target("A+a"),
        )
        .unwrap_err();
        assert_eq!(beyond.reason, InvalidRecordReason::InvalidMm);
    }

    #[test]
    fn parser_accepts_spec_permitted_empty_coordinate_group() {
        let parsed = parse_mod_tags(
            b"T+17802.;A+a.,0;",
            &[255],
            4,
            b"ATAA",
            false,
            MmQuestionMarkPolicy::Unknown,
            &target("A+a"),
        )
        .unwrap();

        assert_eq!(parsed.ml_values_consumed, 1);
        assert_eq!(parsed.calls.len(), 1);
        assert_eq!(parsed.calls[0].query_index, 0);
        assert_eq!(parsed.implicit_calls.len(), 2);
    }

    #[test]
    fn skip_markers_expand_to_distinct_observation_states() {
        let omitted = parse_mod_tags(
            b"A+a,0;",
            &[1],
            2,
            b"AA",
            false,
            MmQuestionMarkPolicy::Unknown,
            &target("A+a"),
        )
        .unwrap();
        assert_eq!(omitted.target_groups_skip_flag_omitted, 1);
        assert_eq!(omitted.implicit_low_probability_candidates, 1);
        assert_eq!(omitted.implicit_calls.len(), 1);
        assert_eq!(
            omitted.implicit_calls[0].state,
            ObservationState::ImplicitBelowEmissionThreshold
        );

        let low = parse_mod_tags(
            b"A+a.,0;",
            &[1],
            2,
            b"AA",
            false,
            MmQuestionMarkPolicy::Unknown,
            &target("A+a"),
        )
        .unwrap();
        assert_eq!(low.target_groups_low_probability, 1);
        assert_eq!(low.implicit_low_probability_candidates, 1);
        assert_eq!(
            low.implicit_calls[0].state,
            ObservationState::ImplicitBelowEmissionThreshold
        );

        let unknown = parse_mod_tags(
            b"A+a?,0;",
            &[1],
            2,
            b"AA",
            false,
            MmQuestionMarkPolicy::Unknown,
            &target("A+a"),
        )
        .unwrap();
        assert_eq!(unknown.target_groups_unknown, 1);
        assert_eq!(unknown.unknown_candidates, 1);
        assert_eq!(unknown.implicit_calls[0].state, ObservationState::Unknown);

        let threshold_sparse = parse_mod_tags(
            b"A+a?,0;",
            &[1],
            2,
            b"AA",
            false,
            MmQuestionMarkPolicy::BelowEmissionThreshold,
            &target("A+a"),
        )
        .unwrap();
        assert_eq!(threshold_sparse.target_groups_unknown, 1);
        assert_eq!(threshold_sparse.unknown_candidates, 0);
        assert_eq!(threshold_sparse.implicit_low_probability_candidates, 1);
        assert_eq!(
            threshold_sparse.implicit_calls[0].state,
            ObservationState::ImplicitBelowEmissionThreshold
        );
    }

    #[test]
    fn normalization_applies_flag_and_mm_strand_and_skips_unmapped_query_bases() {
        let projection = QueryReferenceProjection {
            query_bases: vec![
                QueryBaseProjection::SoftClip,
                QueryBaseProjection::Reference(10),
                QueryBaseProjection::Insertion,
                QueryBaseProjection::Reference(11),
            ],
            reference_end0: 12,
        };
        let parsed = ParsedModTags {
            calls: vec![
                TargetCall {
                    query_index: 0,
                    ml_byte: 1,
                    strand: MmStrand::Forward,
                },
                TargetCall {
                    query_index: 1,
                    ml_byte: 2,
                    strand: MmStrand::Forward,
                },
                TargetCall {
                    query_index: 2,
                    ml_byte: 3,
                    strand: MmStrand::Forward,
                },
            ],
            target_groups_skip_flag_omitted: 1,
            ml_values_consumed: 3,
            ..ParsedModTags::default()
        };
        let imported = normalize_target_calls(
            "chr1",
            "read1",
            true,
            &projection,
            parsed,
            &options(),
            &target("A+a"),
        )
        .unwrap();

        assert_eq!(imported.target_calls_in_soft_clips, 1);
        assert_eq!(imported.target_calls_in_insertions, 1);
        assert_eq!(imported.observations.len(), 1);
        assert_eq!(imported.observations[0].key.site.pos0, 10);
        assert_eq!(imported.observations[0].key.site.strand, Strand::Minus);
        assert_eq!(imported.observations[0].key.read_id, "S1::read1");
        assert_eq!(genomic_strand(true, MmStrand::Reverse), Strand::Plus);
    }

    fn temp_bam(label: &str) -> PathBuf {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        std::env::temp_dir().join(format!(
            "trackcluster-modbam-{label}-{}-{nonce}.bam",
            std::process::id()
        ))
    }

    fn sam_record_with_sequence(
        name: &str,
        flags: Flags,
        mn: usize,
        sequence: &[u8],
    ) -> sam::alignment::RecordBuf {
        use sam::alignment::record_buf::{
            data::field::{value::Array as BufArray, Value as BufValue},
            Cigar, Data, Sequence,
        };

        let cigar: Cigar = [Op::new(Kind::Match, 4)].into_iter().collect();
        let data: Data = [
            (Tag::BASE_MODIFICATIONS, BufValue::from("A+a.,0;")),
            (
                Tag::BASE_MODIFICATION_PROBABILITIES,
                BufValue::Array(BufArray::UInt8(vec![255])),
            ),
            (
                Tag::BASE_MODIFICATION_SEQUENCE_LENGTH,
                BufValue::UInt32(u32::try_from(mn).unwrap()),
            ),
        ]
        .into_iter()
        .collect();

        sam::alignment::RecordBuf::builder()
            .set_name(name)
            .set_flags(flags)
            .set_reference_sequence_id(0)
            .set_alignment_start("101".parse().unwrap())
            .set_mapping_quality(sam::alignment::record::MappingQuality::new(60).unwrap())
            .set_cigar(cigar)
            .set_sequence(Sequence::from(sequence.to_vec()))
            .set_data(data)
            .build()
    }

    fn sam_record(name: &str, flags: Flags, mn: usize) -> sam::alignment::RecordBuf {
        sam_record_with_sequence(name, flags, mn, b"ACCC")
    }

    fn write_test_bam(path: &Path, records: &[sam::alignment::RecordBuf]) {
        use sam::header::record::value::{map::ReferenceSequence, Map};

        let header = sam::Header::builder()
            .add_reference_sequence(
                "chr1",
                Map::<ReferenceSequence>::new(NonZero::new(1000).unwrap()),
            )
            .build();
        let file = File::create(path).unwrap();
        let mut writer = bam::io::Writer::new(file);
        writer.write_header(&header).unwrap();
        for record in records {
            writer.write_alignment_record(&header, record).unwrap();
        }
        writer.try_finish().unwrap();
    }

    #[test]
    fn bam_scan_filters_nonprimary_and_discards_all_duplicate_primaries() {
        let path = temp_bam("duplicates");
        let records = vec![
            sam_record("dup", Flags::empty(), 4),
            sam_record("dup", Flags::empty(), 4),
            sam_record("keep", Flags::empty(), 4),
            sam_record("keep", Flags::SECONDARY, 4),
            sam_record("keep", Flags::SUPPLEMENTARY, 4),
            sam::alignment::RecordBuf::builder()
                .set_name("unmapped")
                .build(),
        ];
        write_test_bam(&path, &records);

        let mut options = options();
        options.invalid_record_policy = InvalidRecordPolicy::Skip;
        options.source_emission_threshold = Some(0.0);
        let result = read_modbam(&path, &options).unwrap();
        assert_eq!(result.qc.total_records, 6);
        assert_eq!(result.qc.skipped_unmapped, 1);
        assert_eq!(result.qc.skipped_secondary, 1);
        assert_eq!(result.qc.skipped_supplementary, 1);
        assert_eq!(result.qc.duplicate_primary_reads, 1);
        assert_eq!(result.qc.skipped_duplicate_primary_records, 2);
        assert_eq!(result.qc.retained_records, 1);
        assert_eq!(result.observations.len(), 1);
        assert_eq!(result.observations[0].key.read_id, "S1::keep");
        assert!(!result.semantics.candidate_observations_complete);

        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn bam_scan_fail_policy_reports_duplicate_primary() {
        let path = temp_bam("duplicate-fail");
        write_test_bam(
            &path,
            &[
                sam_record("dup", Flags::empty(), 4),
                sam_record("dup", Flags::empty(), 4),
            ],
        );

        let error = read_modbam(&path, &options()).unwrap_err();
        assert!(matches!(
            error,
            ModBamError::InvalidRecord {
                reason: InvalidRecordReason::DuplicatePrimary,
                ..
            }
        ));
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn bam_scan_skip_policy_counts_mn_error_and_marks_result_incomplete() {
        let path = temp_bam("invalid-mn");
        write_test_bam(
            &path,
            &[
                sam_record("bad", Flags::empty(), 3),
                sam_record("good", Flags::empty(), 4),
            ],
        );

        let mut options = options();
        options.invalid_record_policy = InvalidRecordPolicy::Skip;
        let result = read_modbam(&path, &options).unwrap();
        assert_eq!(
            result
                .qc
                .invalid_records
                .get(&InvalidRecordReason::MnMismatch),
            Some(&1)
        );
        assert_eq!(result.observations.len(), 1);
        assert!(!result.semantics.candidate_observations_complete);
        assert_eq!(
            result.semantics.implicit_skip_policy,
            ImplicitSkipPolicy::NotApplicable
        );
        assert!(result.semantics.explicit_observations_only);

        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn clean_threshold_zero_input_is_marked_complete() {
        let path = temp_bam("complete");
        write_test_bam(&path, &[sam_record("read1", Flags::empty(), 4)]);

        let mut options = options();
        options.source_emission_threshold = Some(0.0);
        let result = read_modbam(&path, &options).unwrap();
        assert!(result.semantics.candidate_observations_complete);
        assert_eq!(
            result.semantics.implicit_skip_policy,
            ImplicitSkipPolicy::NotApplicable
        );
        assert!(result.qc.invalid_records.is_empty());

        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn implicit_candidates_are_emitted_and_complete_when_projectable() {
        let path = temp_bam("implicit");
        write_test_bam(
            &path,
            &[sam_record_with_sequence(
                "read1",
                Flags::empty(),
                4,
                b"AAAA",
            )],
        );

        let mut options = options();
        options.source_emission_threshold = Some(0.0);
        let result = read_modbam(&path, &options).unwrap();
        assert_eq!(result.qc.implicit_low_probability_candidates, 3);
        assert_eq!(
            result.semantics.implicit_skip_policy,
            ImplicitSkipPolicy::LowProbability
        );
        assert!(result.semantics.candidate_observations_complete);
        assert!(!result.semantics.explicit_observations_only);
        assert_eq!(result.observations.len(), 4);
        assert_eq!(result.qc.emitted_explicit_observations, 1);
        assert_eq!(result.qc.emitted_implicit_observations, 3);

        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn bam_scan_strictly_rejects_non_byte_ml_arrays() {
        use sam::alignment::record_buf::data::field::{
            value::Array as BufArray, Value as BufValue,
        };

        let path = temp_bam("wrong-ml-type");
        let mut record = sam_record("read1", Flags::empty(), 4);
        *record
            .data_mut()
            .get_mut(&Tag::BASE_MODIFICATION_PROBABILITIES)
            .unwrap() = BufValue::Array(BufArray::UInt16(vec![255]));
        write_test_bam(&path, &[record]);

        let mut options = options();
        options.invalid_record_policy = InvalidRecordPolicy::Skip;
        let result = read_modbam(&path, &options).unwrap();
        assert_eq!(
            result
                .qc
                .invalid_records
                .get(&InvalidRecordReason::InvalidMlType),
            Some(&1)
        );
        assert!(result.observations.is_empty());

        let _ = std::fs::remove_file(path);
    }
}
