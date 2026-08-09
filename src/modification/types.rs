use std::fmt;
use std::str::FromStr;

use serde::{Deserialize, Serialize};

use crate::model::Strand;

/// Schema version for normalized modification observations and assay metadata.
pub const MODIFICATION_SCHEMA_VERSION: u32 = 1;

/// Caller interpretation of one candidate read-site observation.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ObservationState {
    /// The caller emitted a numeric read-level probability.
    ExplicitProbability,
    /// The caller guarantees that an omitted candidate is below its emission threshold.
    ImplicitBelowEmissionThreshold,
    /// The candidate exists, but its modification state cannot be interpreted.
    Unknown,
}

impl ObservationState {
    /// Return the stable TSV token.
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::ExplicitProbability => "explicit_probability",
            Self::ImplicitBelowEmissionThreshold => "implicit_below_emission_threshold",
            Self::Unknown => "unknown",
        }
    }
}

impl fmt::Display for ObservationState {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.as_str())
    }
}

impl FromStr for ObservationState {
    type Err = String;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value {
            "explicit_probability" => Ok(Self::ExplicitProbability),
            "implicit_below_emission_threshold" => Ok(Self::ImplicitBelowEmissionThreshold),
            "unknown" => Ok(Self::Unknown),
            _ => Err(format!(
                "invalid observation_state {value:?}; expected explicit_probability, implicit_below_emission_threshold, or unknown"
            )),
        }
    }
}

/// Meaning of candidate positions omitted from a caller's sparse encoding.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ImplicitSkipPolicy {
    /// Omitted candidates are known to be below the source emission threshold.
    LowProbability,
    /// Omitted candidates have unknown state.
    Unknown,
    /// The source format does not use implicit candidates.
    NotApplicable,
}

impl ImplicitSkipPolicy {
    /// Return the stable metadata/QC token.
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::LowProbability => "low_probability",
            Self::Unknown => "unknown",
            Self::NotApplicable => "not_applicable",
        }
    }
}

impl fmt::Display for ImplicitSkipPolicy {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.as_str())
    }
}

/// Structural relationship between an isoform and a genomic modification site.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum SiteState {
    /// The genomic base and required context are present in the isoform.
    Present,
    /// The genomic base is absent from the isoform's exons.
    StructurallyAbsent,
    /// The base is present, but the required sequence context crosses a splice junction.
    ContextDependent,
    /// The strand-oriented reference base does not match the canonical modified base.
    ReferenceBaseMismatch,
    /// A source transcript coordinate could not be projected reliably.
    Unprojectable,
}

impl SiteState {
    /// Return the stable output token.
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Present => "present",
            Self::StructurallyAbsent => "structurally_absent",
            Self::ContextDependent => "context_dependent",
            Self::ReferenceBaseMismatch => "reference_base_mismatch",
            Self::Unprojectable => "unprojectable",
        }
    }
}

impl fmt::Display for SiteState {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.as_str())
    }
}

/// Evidence used to count reads covering a genomic base.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum CoverageBasis {
    /// No independent read coverage source was available.
    Unavailable,
    /// Coverage was approximated from BED exon spans.
    BedApproximate,
    /// Coverage was calculated from base-level BAM/CIGAR projection.
    BamExact,
}

impl CoverageBasis {
    /// Return the stable output token.
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Unavailable => "unavailable",
            Self::BedApproximate => "bed_approximate",
            Self::BamExact => "bam_exact",
        }
    }
}

impl fmt::Display for CoverageBasis {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.as_str())
    }
}

/// Stable reason explaining whether a site-level row is analysis eligible.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum EligibilityReason {
    /// Every V1 eligibility condition is satisfied.
    Ok,
    /// Too few callable molecules passed the configured minimum.
    LowCallable,
    /// The site is structurally absent from this isoform.
    SiteAbsent,
    /// The candidate context depends on splice structure.
    ContextDependent,
    /// The strand-oriented reference base does not match the canonical modified base.
    ReferenceBaseMismatch,
    /// The caller output does not contain a complete retained candidate universe.
    IncompleteCandidateUniverse,
    /// Observation-to-assignment join rate is below the configured minimum.
    JoinRateLow,
    /// Observation-to-assignment join rate at this genomic site is below the minimum.
    SiteJoinRateLow,
    /// At least one candidate state needed for the fraction denominator is unknown.
    UnknownDenominator,
    /// The site could not be projected to a genomic base.
    Unprojectable,
}

impl EligibilityReason {
    /// Return the stable output token.
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Ok => "ok",
            Self::LowCallable => "low_callable",
            Self::SiteAbsent => "site_absent",
            Self::ContextDependent => "context_dependent",
            Self::ReferenceBaseMismatch => "reference_base_mismatch",
            Self::IncompleteCandidateUniverse => "incomplete_candidate_universe",
            Self::JoinRateLow => "join_rate_low",
            Self::SiteJoinRateLow => "site_join_rate_low",
            Self::UnknownDenominator => "unknown_denominator",
            Self::Unprojectable => "unprojectable",
        }
    }

    /// Return whether the reason denotes an eligible row.
    pub const fn is_eligible(self) -> bool {
        matches!(self, Self::Ok)
    }
}

impl fmt::Display for EligibilityReason {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.as_str())
    }
}

/// Genomic site identity shared across isoforms.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ModSiteKey {
    /// Reference sequence name.
    pub chrom: String,
    /// Zero-based genomic base coordinate.
    pub pos0: u32,
    /// Genomic strand.
    pub strand: Strand,
    /// SAM-style canonical base, strand, and modification code (for example `A+a`).
    pub mod_code: String,
}

impl ModSiteKey {
    /// Return the stable genomic site identifier, excluding modification code.
    pub fn site_id(&self) -> String {
        format!("{}:{}:{}", self.chrom, self.pos0, self.strand.as_char())
    }
}

/// Uniqueness key for one normalized read-site observation.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ModObservationKey {
    /// Assay compatibility stratum.
    pub assay_id: String,
    /// Biological sample identifier.
    pub sample: String,
    /// Final TrackCluster molecule identifier.
    pub read_id: String,
    /// Genomic site and modification identity.
    pub site: ModSiteKey,
}

/// One model-candidate read-site observation.
#[derive(Clone, Debug, PartialEq)]
pub struct ModObservation {
    /// Assay/sample/read/site identity.
    pub key: ModObservationKey,
    /// Numeric caller probability for explicit observations.
    pub probability: Option<f64>,
    /// Caller observation semantics.
    pub observation_state: ObservationState,
    /// Caller sequence context or candidate-rule label.
    pub context: Option<String>,
    /// Source transcript identifier retained for provenance.
    pub source_transcript_id: Option<String>,
    /// Zero-based source transcript offset retained for provenance.
    pub source_pos0: Option<u64>,
}

impl ModObservation {
    /// Validate the cross-field normalized observation contract.
    pub fn validate(&self) -> Result<(), String> {
        validate_identifier("assay_id", &self.key.assay_id)?;
        validate_identifier("sample", &self.key.sample)?;
        validate_identifier("read_id", &self.key.read_id)?;
        validate_identifier("chrom", &self.key.site.chrom)?;
        validate_mod_code(&self.key.site.mod_code)?;

        let read_sample = crate::sample::split_tagged_read_name(&self.key.read_id)
            .ok_or_else(|| {
                format!(
                    "read_id {:?} must use the <sample>::<read_id> form",
                    self.key.read_id
                )
            })?
            .0;
        if read_sample != self.key.sample {
            return Err(format!(
                "read_id {:?} has sample prefix {read_sample:?}, expected {:?}",
                self.key.read_id, self.key.sample
            ));
        }

        match (self.observation_state, self.probability) {
            (ObservationState::ExplicitProbability, Some(probability))
                if probability.is_finite() && (0.0..=1.0).contains(&probability) => {}
            (ObservationState::ExplicitProbability, Some(probability)) => {
                return Err(format!(
                    "explicit probability must be finite and in [0, 1], got {probability}"
                ));
            }
            (ObservationState::ExplicitProbability, None) => {
                return Err("explicit_probability requires a probability".to_owned());
            }
            (_, Some(_)) => {
                return Err(format!(
                    "{} requires probability=NA",
                    self.observation_state
                ));
            }
            (_, None) => {}
        }

        if let Some(context) = &self.context {
            validate_identifier("context", context)?;
        }
        if let Some(transcript_id) = &self.source_transcript_id {
            validate_identifier("source_transcript_id", transcript_id)?;
        }
        if self.source_transcript_id.is_some() != self.source_pos0.is_some() {
            return Err(
                "source_transcript_id and source_pos0 must either both be present or both be NA"
                    .to_owned(),
            );
        }
        Ok(())
    }
}

/// Dataset-level provenance for normalized modification observations.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AssayMetadata {
    /// Metadata schema version; V1 requires `1`.
    pub schema_version: u32,
    /// Compatibility stratum used by aggregation and contrasts.
    pub assay_id: String,
    /// Caller name, such as `dorado` or `m6anet`.
    pub caller: String,
    /// Caller version, or `unknown` when it cannot be recovered.
    pub caller_version: String,
    /// Exact model identifier.
    pub model_id: String,
    /// Sequencing chemistry, such as `RNA002` or `RNA004`.
    pub chemistry: String,
    /// Candidate universe rule, such as `DRACH` or `all-context-A`.
    pub candidate_rule: String,
    /// Probability threshold below which the source omitted candidates.
    pub source_emission_threshold: Option<f64>,
    /// Caller/site filtering provenance.
    pub source_site_filter: String,
    /// Whether all observations in the retained candidate universe are represented.
    pub candidate_observations_complete: bool,
    /// Interpretation of candidates omitted by the source encoding.
    pub implicit_skip_policy: ImplicitSkipPolicy,
    /// How source coordinates became genomic coordinates.
    pub coordinate_source: String,
    /// How source read identities became final TrackCluster read IDs.
    pub read_id_mapping: String,
    /// Source paths or immutable source identifiers.
    #[serde(default)]
    pub source_files: Vec<String>,
}

impl AssayMetadata {
    /// Validate metadata required for reproducible threshold interpretation.
    pub fn validate(&self) -> Result<(), String> {
        if self.schema_version != MODIFICATION_SCHEMA_VERSION {
            return Err(format!(
                "unsupported assay metadata schema_version {}; expected {}",
                self.schema_version, MODIFICATION_SCHEMA_VERSION
            ));
        }
        for (field, value) in [
            ("assay_id", self.assay_id.as_str()),
            ("caller", self.caller.as_str()),
            ("caller_version", self.caller_version.as_str()),
            ("model_id", self.model_id.as_str()),
            ("chemistry", self.chemistry.as_str()),
            ("candidate_rule", self.candidate_rule.as_str()),
            ("source_site_filter", self.source_site_filter.as_str()),
            ("coordinate_source", self.coordinate_source.as_str()),
            ("read_id_mapping", self.read_id_mapping.as_str()),
        ] {
            validate_identifier(field, value)?;
        }
        if let Some(threshold) = self.source_emission_threshold {
            if !threshold.is_finite() || !(0.0..=1.0).contains(&threshold) {
                return Err(format!(
                    "source_emission_threshold must be finite and in [0, 1], got {threshold}"
                ));
            }
        }
        if self.implicit_skip_policy == ImplicitSkipPolicy::LowProbability
            && self.source_emission_threshold.is_none()
        {
            return Err(
                "implicit_skip_policy=low_probability requires source_emission_threshold"
                    .to_owned(),
            );
        }
        for source in &self.source_files {
            validate_identifier("source_files entry", source)?;
        }
        Ok(())
    }
}

pub(crate) fn validate_identifier(field: &str, value: &str) -> Result<(), String> {
    if value.trim().is_empty() || value == "NA" {
        return Err(format!("{field} must not be empty or NA"));
    }
    if value.chars().any(char::is_control) {
        return Err(format!("{field} contains a control character"));
    }
    Ok(())
}

pub(crate) fn iupac_dna_mask(base: u8) -> Option<u8> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0b0001),
        b'C' => Some(0b0010),
        b'G' => Some(0b0100),
        b'T' | b'U' => Some(0b1000),
        b'R' => Some(0b0101),
        b'Y' => Some(0b1010),
        b'S' => Some(0b0110),
        b'W' => Some(0b1001),
        b'K' => Some(0b1100),
        b'M' => Some(0b0011),
        b'B' => Some(0b1110),
        b'D' => Some(0b1101),
        b'H' => Some(0b1011),
        b'V' => Some(0b0111),
        b'N' => Some(0b1111),
        _ => None,
    }
}

pub(crate) fn normalize_centered_iupac_motif(
    value: &str,
    target_base: u8,
) -> Result<String, String> {
    if value.len() < 3 || value.len().is_multiple_of(2) || value.len() > 31 {
        return Err(format!(
            "candidate_rule {value:?} must be an odd-length centered IUPAC DNA motif of 3..31 bases"
        ));
    }

    let mut normalized = String::with_capacity(value.len());
    for raw_base in value.bytes() {
        let base = match raw_base.to_ascii_uppercase() {
            b'U' => b'T',
            base => base,
        };
        if iupac_dna_mask(base).is_none() {
            return Err(format!(
                "candidate_rule {value:?} contains unsupported IUPAC base {:?}",
                char::from(raw_base)
            ));
        }
        normalized.push(char::from(base));
    }

    let target_base = match target_base.to_ascii_uppercase() {
        b'U' => b'T',
        base => base,
    };
    let center = normalized.len() / 2;
    if normalized.as_bytes()[center] != target_base {
        return Err(format!(
            "candidate_rule {value:?} must have target canonical base {} at its center",
            char::from(target_base)
        ));
    }

    Ok(normalized)
}

pub(crate) fn sequence_matches_iupac_motif(sequence: &str, motif: &str) -> bool {
    sequence.len() == motif.len()
        && sequence
            .bytes()
            .zip(motif.bytes())
            .all(|(base, motif_base)| {
                let sequence_mask = match base.to_ascii_uppercase() {
                    b'A' => Some(0b0001),
                    b'C' => Some(0b0010),
                    b'G' => Some(0b0100),
                    b'T' | b'U' => Some(0b1000),
                    _ => None,
                };
                sequence_mask.is_some_and(|mask| {
                    iupac_dna_mask(motif_base).is_some_and(|motif_mask| mask & motif_mask != 0)
                })
            })
}

fn validate_mod_code(value: &str) -> Result<(), String> {
    validate_identifier("mod_code", value)?;
    let bytes = value.as_bytes();
    if bytes.len() < 3
        || !matches!(bytes[0], b'A' | b'C' | b'G' | b'T' | b'U' | b'N')
        || !matches!(bytes[1], b'+' | b'-')
        || !bytes[2..].iter().all(u8::is_ascii_alphanumeric)
    {
        return Err(format!(
            "invalid mod_code {value:?}; expected SAM-style form such as A+a or C+76792"
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn observation(state: ObservationState, probability: Option<f64>) -> ModObservation {
        ModObservation {
            key: ModObservationKey {
                assay_id: "rna004_m6a".to_owned(),
                sample: "S1".to_owned(),
                read_id: "S1::read-1".to_owned(),
                site: ModSiteKey {
                    chrom: "chr1".to_owned(),
                    pos0: 10,
                    strand: Strand::Plus,
                    mod_code: "A+a".to_owned(),
                },
            },
            probability,
            observation_state: state,
            context: Some("DRACH".to_owned()),
            source_transcript_id: None,
            source_pos0: None,
        }
    }

    #[test]
    fn observation_state_controls_probability_presence() {
        assert!(
            observation(ObservationState::ExplicitProbability, Some(0.0))
                .validate()
                .is_ok()
        );
        assert!(
            observation(ObservationState::ExplicitProbability, Some(1.0))
                .validate()
                .is_ok()
        );
        assert!(observation(ObservationState::ExplicitProbability, None)
            .validate()
            .is_err());
        assert!(observation(ObservationState::Unknown, Some(0.0))
            .validate()
            .is_err());
        assert!(observation(ObservationState::Unknown, None)
            .validate()
            .is_ok());
    }

    #[test]
    fn observation_rejects_invalid_probability_and_sample_prefix() {
        assert!(
            observation(ObservationState::ExplicitProbability, Some(f64::NAN))
                .validate()
                .is_err()
        );
        assert!(
            observation(ObservationState::ExplicitProbability, Some(1.1))
                .validate()
                .is_err()
        );
        let mut value = observation(ObservationState::ExplicitProbability, Some(0.5));
        value.key.read_id = "S2::read-1".to_owned();
        assert!(value.validate().is_err());
    }

    #[test]
    fn low_probability_skip_policy_requires_threshold() {
        let metadata = AssayMetadata {
            schema_version: 1,
            assay_id: "rna004_m6a".to_owned(),
            caller: "dorado".to_owned(),
            caller_version: "2.0.0".to_owned(),
            model_id: "rna004_hac_m6a".to_owned(),
            chemistry: "RNA004".to_owned(),
            candidate_rule: "all-context-A".to_owned(),
            source_emission_threshold: None,
            source_site_filter: "none".to_owned(),
            candidate_observations_complete: true,
            implicit_skip_policy: ImplicitSkipPolicy::LowProbability,
            coordinate_source: "genome_aligned_bam".to_owned(),
            read_id_mapping: "bam_qname_with_sample_prefix".to_owned(),
            source_files: Vec::new(),
        };

        assert!(metadata.validate().is_err());
        assert!(AssayMetadata {
            source_emission_threshold: Some(0.0),
            ..metadata
        }
        .validate()
        .is_ok());
    }
}
