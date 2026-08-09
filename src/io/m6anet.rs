//! Strict m6Anet RNA002 read-level probability import.

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::io::Read;
use std::path::Path;

use anyhow::Context;
use csv::{Reader, ReaderBuilder, StringRecord};

use crate::io::gff::{open_maybe_gzip, read_annotation_transcripts_where, GffToBiggOptions};
use crate::model::{Strand, Transcript};
use crate::modification::types::MODIFICATION_SCHEMA_VERSION;
use crate::modification::{
    AssayMetadata, ImplicitSkipPolicy, ModObservation, ModObservationKey, ModSiteKey,
    ObservationState,
};
use crate::sample::{split_tagged_read_name, tagged_read_name, SAMPLE_DELIM};

const INDIVIDUAL_COLUMNS: [&str; 4] = [
    "transcript_id",
    "transcript_position",
    "read_index",
    "probability_modified",
];
const READ_MAP_COLUMNS: [&str; 2] = ["read_index", "read_id"];
const DATA_INFO_COLUMNS: [&str; 5] = [
    "transcript_id",
    "transcript_position",
    "start",
    "end",
    "n_reads",
];
const SITE_PROBABILITY_COLUMNS: [&str; 6] = [
    "transcript_id",
    "transcript_position",
    "n_reads",
    "probability_modified",
    "kmer",
    "mod_ratio",
];
const M6A_MOD_CODE: &str = "A+a";
const SITE_PROBABILITY_THRESHOLD: f64 = 0.9;

/// Official m6Anet read-probability threshold for the HCT116 RNA002 model.
pub const HCT116_RNA002_READ_THRESHOLD: f64 = 0.033_379_376;
/// Official m6Anet read-probability threshold for the Arabidopsis RNA002 model.
pub const ARABIDOPSIS_RNA002_READ_THRESHOLD: f64 = 0.003_297_804_621_979_6;

/// Return the official read-level threshold preset for a known RNA002 model.
pub fn read_probability_threshold_preset(model_id: &str) -> Option<f64> {
    match model_id.to_ascii_lowercase().as_str() {
        "hct116_rna002" => Some(HCT116_RNA002_READ_THRESHOLD),
        "arabidopsis" | "arabidopsis_rna002" => Some(ARABIDOPSIS_RNA002_READ_THRESHOLD),
        _ => None,
    }
}

/// Configuration and provenance for one m6Anet RNA002 import.
#[derive(Clone, Debug, PartialEq)]
pub struct M6anetImportOptions {
    /// Biological sample identifier used in normalized observations.
    pub sample: String,
    /// Assay compatibility stratum used in normalized observations.
    pub assay_id: String,
    /// m6Anet version, or `unknown` when unavailable.
    pub caller_version: String,
    /// Exact m6Anet model identifier.
    pub model_id: String,
    /// Candidate motif rule represented by the source model.
    pub candidate_rule: String,
    /// GTF/GFF parsing options.
    pub annotation_options: GffToBiggOptions,
    /// Optional m6Anet minimum-read threshold used to audit `data.info`.
    pub minimum_reads: Option<u64>,
    /// Source read threshold used only to cross-check `data.site_proba` mod ratios.
    pub read_probability_threshold: Option<f64>,
}

impl M6anetImportOptions {
    /// Construct RNA002 options with strict defaults and explicit sample/assay/model IDs.
    pub fn new(
        sample: impl Into<String>,
        assay_id: impl Into<String>,
        model_id: impl Into<String>,
    ) -> Self {
        let model_id = model_id.into();
        Self {
            sample: sample.into(),
            assay_id: assay_id.into(),
            caller_version: "unknown".to_owned(),
            read_probability_threshold: read_probability_threshold_preset(&model_id),
            model_id,
            candidate_rule: "DRACH".to_owned(),
            annotation_options: GffToBiggOptions::default(),
            minimum_reads: None,
        }
    }
}

/// Audit details derived from an optional m6Anet `data.info` file.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct M6anetDataInfoQc {
    /// Unique transcript sites listed by `data.info`.
    pub sites: usize,
    /// Sites represented by `data.indiv_proba.csv`.
    pub retained_sites: usize,
    /// `data.info` sites absent from `data.indiv_proba.csv`.
    pub filtered_sites: usize,
    /// Sum of `n_reads` over all `data.info` sites.
    pub total_reads: u64,
    /// Sum of `n_reads` over retained sites.
    pub retained_reads: u64,
    /// Minimum-read threshold used for exact retained-site auditing.
    pub minimum_reads: Option<u64>,
}

/// Cross-check counters from an optional m6Anet `data.site_proba.csv` file.
#[derive(Clone, Debug, PartialEq)]
pub struct M6anetSiteProbabilityQc {
    /// Unique retained sites cross-checked against individual probabilities.
    pub sites: usize,
    /// Total read observations represented by those sites.
    pub total_reads: u64,
    /// Sites meeting m6Anet's separate recommended site-probability threshold.
    pub sites_at_or_above_probability_threshold: usize,
    /// Site-level probability threshold used only for this QC count.
    pub site_probability_threshold: f64,
    /// Read threshold used to reproduce `mod_ratio`, when known.
    pub read_probability_threshold: Option<f64>,
}

/// Deterministic counters produced by an m6Anet import.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct M6anetImportQc {
    /// Physical rows read from `data.indiv_proba.csv`.
    pub input_rows: usize,
    /// Unique normalized observations returned.
    pub unique_observations: usize,
    /// Exact duplicate observations folded during import.
    pub duplicate_exact: usize,
    /// Entries in the explicit read-index map.
    pub read_map_entries: usize,
    /// Read-index map entries referenced by the input probabilities.
    pub read_map_entries_used: usize,
    /// Distinct source transcript IDs represented by retained observations.
    pub source_transcripts: usize,
    /// Annotation transcripts materialized for source-coordinate projection.
    pub projection_transcripts_loaded: usize,
    /// Distinct source transcript-position pairs represented by retained observations.
    pub source_sites: usize,
    /// Optional `data.info` audit result.
    pub data_info: Option<M6anetDataInfoQc>,
    /// Optional `data.site_proba.csv` cross-check result.
    pub site_probability: Option<M6anetSiteProbabilityQc>,
}

/// Normalized observations, assay provenance, and import audit counters.
#[derive(Clone, Debug, PartialEq)]
pub struct M6anetImportResult {
    /// Unique observations in deterministic key order.
    pub observations: Vec<ModObservation>,
    /// Assay metadata suitable for the normalized modification sidecar.
    pub metadata: AssayMetadata,
    /// Import counters and optional `data.info` audit.
    pub qc: M6anetImportQc,
}

#[derive(Clone, Debug)]
struct ReadMapEntry {
    read_id: String,
    line: usize,
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct TranscriptSite {
    transcript_id: String,
    pos0: u64,
}

#[derive(Clone, Copy, Debug)]
struct DataInfoEntry {
    line: usize,
    n_reads: u64,
}

#[derive(Debug)]
struct TranscriptIndex {
    exact: HashMap<String, usize>,
    versions: BTreeMap<String, Vec<String>>,
}

impl TranscriptIndex {
    fn new(transcripts: &[Transcript]) -> anyhow::Result<Self> {
        let mut exact = HashMap::with_capacity(transcripts.len());
        let mut versions: BTreeMap<String, Vec<String>> = BTreeMap::new();
        for (index, transcript) in transcripts.iter().enumerate() {
            if exact.insert(transcript.name.clone(), index).is_some() {
                anyhow::bail!(
                    "annotation contains duplicate transcript ID {:?}",
                    transcript.name
                );
            }
            versions
                .entry(unversioned_transcript_id(&transcript.name).to_owned())
                .or_default()
                .push(transcript.name.clone());
        }
        for values in versions.values_mut() {
            values.sort();
        }
        Ok(Self { exact, versions })
    }

    fn resolve<'a>(
        &self,
        transcripts: &'a [Transcript],
        transcript_id: &str,
    ) -> anyhow::Result<&'a Transcript> {
        if let Some(index) = self.exact.get(transcript_id) {
            return Ok(&transcripts[*index]);
        }

        let base = unversioned_transcript_id(transcript_id);
        if let Some(available) = self.versions.get(base) {
            anyhow::bail!(
                "m6Anet transcript version mismatch for {transcript_id:?}: exact ID is absent; annotation contains {available:?}"
            );
        }
        anyhow::bail!("m6Anet transcript {transcript_id:?} is absent from the annotation")
    }
}

fn unversioned_transcript_id(transcript_id: &str) -> &str {
    match transcript_id.rsplit_once('.') {
        Some((base, version))
            if !base.is_empty()
                && !version.is_empty()
                && version.bytes().all(|byte| byte.is_ascii_digit()) =>
        {
            base
        }
        _ => transcript_id,
    }
}

fn validate_token(field: &str, value: &str) -> anyhow::Result<()> {
    if value.is_empty() || value == "NA" {
        anyhow::bail!("{field} must not be empty or NA");
    }
    if value.trim() != value {
        anyhow::bail!("{field} {value:?} must not have leading or trailing whitespace");
    }
    if value.chars().any(char::is_control) {
        anyhow::bail!("{field} {value:?} contains a control character");
    }
    Ok(())
}

fn validate_options(options: &M6anetImportOptions) -> anyhow::Result<String> {
    for (field, value) in [
        ("sample", options.sample.as_str()),
        ("assay_id", options.assay_id.as_str()),
        ("caller_version", options.caller_version.as_str()),
        ("model_id", options.model_id.as_str()),
        ("candidate_rule", options.candidate_rule.as_str()),
    ] {
        validate_token(field, value)?;
    }
    if options.sample.contains(SAMPLE_DELIM) {
        anyhow::bail!(
            "sample {:?} must not contain the reserved delimiter {SAMPLE_DELIM:?}",
            options.sample
        );
    }
    if options
        .read_probability_threshold
        .is_some_and(|value| !value.is_finite() || !(0.0..=1.0).contains(&value))
    {
        anyhow::bail!("read_probability_threshold must be finite and in [0, 1]");
    }
    crate::modification::types::normalize_centered_iupac_motif(&options.candidate_rule, b'A')
        .map_err(anyhow::Error::msg)
}

fn exact_header<R: Read>(
    reader: &mut Reader<R>,
    expected: &[&str],
    description: &str,
    path: &Path,
) -> anyhow::Result<()> {
    let actual = reader
        .headers()
        .with_context(|| format!("read {description} header {path:?}"))?
        .clone();
    let expected_record = StringRecord::from(expected.to_vec());
    if actual != expected_record {
        anyhow::bail!(
            "{description} {path:?} header mismatch: expected {expected:?}, found {:?}",
            actual.iter().collect::<Vec<_>>()
        );
    }
    Ok(())
}

fn required<'a>(record: &'a StringRecord, index: usize, field: &str) -> anyhow::Result<&'a str> {
    let value = record.get(index).unwrap_or_default();
    if value.is_empty() {
        anyhow::bail!("{field} must not be empty");
    }
    Ok(value)
}

fn parse_u64(record: &StringRecord, index: usize, field: &str) -> anyhow::Result<u64> {
    let value = required(record, index, field)?;
    value
        .parse::<u64>()
        .with_context(|| format!("invalid unsigned integer for {field}: {value:?}"))
}

fn parse_probability(record: &StringRecord, index: usize) -> anyhow::Result<f64> {
    parse_unit_interval(record, index, "probability_modified")
}

fn parse_unit_interval(record: &StringRecord, index: usize, field: &str) -> anyhow::Result<f64> {
    let value = required(record, index, field)?;
    let probability = value
        .parse::<f64>()
        .with_context(|| format!("invalid {field} {value:?}"))?;
    if !probability.is_finite() || !(0.0..=1.0).contains(&probability) {
        anyhow::bail!("{field} must be finite and in [0, 1], got {value:?}");
    }
    Ok(probability)
}

fn normalize_read_id(sample: &str, mapped_id: &str) -> anyhow::Result<String> {
    validate_token("read_id", mapped_id)?;
    if mapped_id.contains(SAMPLE_DELIM) {
        let (mapped_sample, _) = split_tagged_read_name(mapped_id).with_context(|| {
            format!("mapped read_id {mapped_id:?} has malformed {SAMPLE_DELIM:?} sample prefix")
        })?;
        if mapped_sample != sample {
            anyhow::bail!(
                "mapped read_id {mapped_id:?} has sample prefix {mapped_sample:?}, expected {sample:?}"
            );
        }
        Ok(mapped_id.to_owned())
    } else {
        Ok(tagged_read_name(sample, mapped_id))
    }
}

fn read_read_map(path: &Path, sample: &str) -> anyhow::Result<BTreeMap<String, ReadMapEntry>> {
    let input = open_maybe_gzip(path).with_context(|| format!("open m6Anet read map {path:?}"))?;
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .flexible(false)
        .from_reader(input);
    exact_header(&mut reader, &READ_MAP_COLUMNS, "m6Anet read map", path)?;

    let mut by_index: BTreeMap<String, ReadMapEntry> = BTreeMap::new();
    let mut by_read_id: HashMap<String, (String, usize)> = HashMap::new();
    for (record_index, result) in reader.records().enumerate() {
        let line = record_index + 2;
        let record = result.with_context(|| format!("parse m6Anet read map {path:?}:{line}"))?;
        if record.len() != READ_MAP_COLUMNS.len() {
            anyhow::bail!(
                "m6Anet read map {path:?}:{line}: expected {} fields, found {}",
                READ_MAP_COLUMNS.len(),
                record.len()
            );
        }
        let read_index = required(&record, 0, "read_index")?;
        validate_token("read_index", read_index)
            .with_context(|| format!("validate m6Anet read map {path:?}:{line}"))?;
        let read_id = normalize_read_id(sample, required(&record, 1, "read_id")?)
            .with_context(|| format!("validate m6Anet read map {path:?}:{line}"))?;
        if let Some(previous) = by_index.get(read_index) {
            anyhow::bail!(
                "m6Anet read map {path:?}:{line}: duplicate read_index {read_index}, first seen at line {}",
                previous.line
            );
        }
        if let Some((previous_index, previous_line)) = by_read_id.get(&read_id) {
            anyhow::bail!(
                "m6Anet read map {path:?}:{line}: read_id {read_id:?} maps from both read_index {previous_index} at line {previous_line} and {read_index}"
            );
        }
        by_read_id.insert(read_id.clone(), (read_index.to_owned(), line));
        by_index.insert(read_index.to_owned(), ReadMapEntry { read_id, line });
    }
    if by_index.is_empty() {
        anyhow::bail!("m6Anet read map {path:?} contains no mappings");
    }
    Ok(by_index)
}

fn read_source_transcript_ids(
    path: &Path,
    expected_columns: &[&str],
    description: &str,
) -> anyhow::Result<BTreeSet<String>> {
    let input = open_maybe_gzip(path).with_context(|| format!("open {description} {path:?}"))?;
    let mut reader = ReaderBuilder::new()
        .delimiter(b',')
        .has_headers(true)
        .flexible(false)
        .from_reader(input);
    exact_header(&mut reader, expected_columns, description, path)?;

    let mut transcript_ids = BTreeSet::new();
    for (record_index, result) in reader.records().enumerate() {
        let line = record_index + 2;
        let record = result.with_context(|| format!("parse {description} {path:?}:{line}"))?;
        if record.len() != expected_columns.len() {
            anyhow::bail!(
                "{description} {path:?}:{line}: expected {} fields, found {}",
                expected_columns.len(),
                record.len()
            );
        }
        let transcript_id = required(&record, 0, "transcript_id")?;
        validate_token("transcript_id", transcript_id)
            .with_context(|| format!("validate {description} {path:?}:{line}"))?;
        transcript_ids.insert(transcript_id.to_owned());
    }
    if transcript_ids.is_empty() {
        anyhow::bail!("{description} {path:?} contains no transcript IDs");
    }
    Ok(transcript_ids)
}

fn project_transcript_position(
    transcript: &Transcript,
    transcript_id: &str,
    pos0: u64,
) -> anyhow::Result<u32> {
    let geometry = transcript.geometry();
    if geometry.strand == Strand::Unknown {
        anyhow::bail!(
            "m6Anet transcript {transcript_id:?} has unknown strand; transcript-oriented projection is undefined"
        );
    }
    let spliced_len = geometry.spliced_len();
    if pos0 >= spliced_len {
        anyhow::bail!(
            "m6Anet transcript_position {pos0} is out of bounds for transcript {transcript_id:?} with spliced length {spliced_len}"
        );
    }
    geometry
        .spliced_offset_to_genomic(pos0)
        .map(|position| position.get())
        .with_context(|| {
            format!("project m6Anet transcript_position {pos0} for transcript {transcript_id:?}")
        })
}

fn audit_data_info(
    path: &Path,
    observed_site_reads: &BTreeMap<TranscriptSite, u64>,
    minimum_reads: Option<u64>,
    transcripts: &[Transcript],
    transcript_index: &TranscriptIndex,
) -> anyhow::Result<M6anetDataInfoQc> {
    let input = open_maybe_gzip(path).with_context(|| format!("open m6Anet data.info {path:?}"))?;
    let mut reader = ReaderBuilder::new()
        .delimiter(b',')
        .has_headers(true)
        .flexible(false)
        .from_reader(input);
    exact_header(&mut reader, &DATA_INFO_COLUMNS, "m6Anet data.info", path)?;

    let mut entries = BTreeMap::new();
    for (record_index, result) in reader.records().enumerate() {
        let line = record_index + 2;
        let record = result.with_context(|| format!("parse m6Anet data.info {path:?}:{line}"))?;
        if record.len() != DATA_INFO_COLUMNS.len() {
            anyhow::bail!(
                "m6Anet data.info {path:?}:{line}: expected {} fields, found {}",
                DATA_INFO_COLUMNS.len(),
                record.len()
            );
        }
        let transcript_id = required(&record, 0, "transcript_id")?;
        validate_token("transcript_id", transcript_id)
            .with_context(|| format!("validate m6Anet data.info {path:?}:{line}"))?;
        let pos0 = parse_u64(&record, 1, "transcript_position")
            .with_context(|| format!("validate m6Anet data.info {path:?}:{line}"))?;
        let start = parse_u64(&record, 2, "start")
            .with_context(|| format!("validate m6Anet data.info {path:?}:{line}"))?;
        let end = parse_u64(&record, 3, "end")
            .with_context(|| format!("validate m6Anet data.info {path:?}:{line}"))?;
        let n_reads = parse_u64(&record, 4, "n_reads")
            .with_context(|| format!("validate m6Anet data.info {path:?}:{line}"))?;
        if end <= start {
            anyhow::bail!(
                "m6Anet data.info {path:?}:{line}: end {end} must be greater than start {start}"
            );
        }
        if n_reads == 0 {
            anyhow::bail!("m6Anet data.info {path:?}:{line}: n_reads must be greater than zero");
        }

        let transcript = transcript_index
            .resolve(transcripts, transcript_id)
            .with_context(|| format!("validate m6Anet data.info {path:?}:{line}"))?;
        project_transcript_position(transcript, transcript_id, pos0)
            .with_context(|| format!("validate m6Anet data.info {path:?}:{line}"))?;

        let site = TranscriptSite {
            transcript_id: transcript_id.to_owned(),
            pos0,
        };
        if let Some(previous) = entries.insert(site.clone(), DataInfoEntry { line, n_reads }) {
            anyhow::bail!(
                "m6Anet data.info {path:?}:{line}: duplicate site {:?}:{} first seen at line {}",
                site.transcript_id,
                site.pos0,
                previous.line
            );
        }
    }
    if entries.is_empty() {
        anyhow::bail!("m6Anet data.info {path:?} contains no sites");
    }

    for (site, observed_reads) in observed_site_reads {
        let entry = entries.get(site).with_context(|| {
            format!(
                "m6Anet data.info {path:?} is missing retained site {:?}:{}",
                site.transcript_id, site.pos0
            )
        })?;
        if *observed_reads != entry.n_reads {
            anyhow::bail!(
                "m6Anet data.info {path:?}: retained site {:?}:{} has n_reads={}, but data.indiv_proba contains {} rows",
                site.transcript_id,
                site.pos0,
                entry.n_reads,
                observed_reads
            );
        }
    }

    if let Some(minimum_reads) = minimum_reads {
        for (site, entry) in &entries {
            let expected_retained = entry.n_reads >= minimum_reads;
            let observed_retained = observed_site_reads.contains_key(site);
            if expected_retained != observed_retained {
                anyhow::bail!(
                    "m6Anet data.info {path:?}: retained-site mismatch for {:?}:{}: n_reads={} with minimum_reads={minimum_reads} implies retained={expected_retained}, but data.indiv_proba retained={observed_retained}",
                    site.transcript_id,
                    site.pos0,
                    entry.n_reads
                );
            }
        }
    }

    let total_reads = entries.values().try_fold(0u64, |sum, entry| {
        sum.checked_add(entry.n_reads)
            .context("m6Anet data.info total n_reads overflow")
    })?;
    let retained_reads = observed_site_reads.keys().try_fold(0u64, |sum, site| {
        sum.checked_add(entries[site].n_reads)
            .context("m6Anet data.info retained n_reads overflow")
    })?;
    Ok(M6anetDataInfoQc {
        sites: entries.len(),
        retained_sites: observed_site_reads.len(),
        filtered_sites: entries.len() - observed_site_reads.len(),
        total_reads,
        retained_reads,
        minimum_reads,
    })
}

fn audit_site_probability(
    path: &Path,
    observed_site_probabilities: &BTreeMap<TranscriptSite, Vec<f64>>,
    read_probability_threshold: Option<f64>,
    candidate_rule: &str,
) -> anyhow::Result<M6anetSiteProbabilityQc> {
    let input =
        open_maybe_gzip(path).with_context(|| format!("open m6Anet data.site_proba {path:?}"))?;
    let mut reader = ReaderBuilder::new()
        .delimiter(b',')
        .has_headers(true)
        .flexible(false)
        .from_reader(input);
    exact_header(
        &mut reader,
        &SITE_PROBABILITY_COLUMNS,
        "m6Anet data.site_proba",
        path,
    )?;

    let mut seen = BTreeSet::new();
    let mut total_reads = 0u64;
    let mut sites_at_or_above_probability_threshold = 0usize;
    for (record_index, result) in reader.records().enumerate() {
        let line = record_index + 2;
        let record =
            result.with_context(|| format!("parse m6Anet data.site_proba {path:?}:{line}"))?;
        if record.len() != SITE_PROBABILITY_COLUMNS.len() {
            anyhow::bail!(
                "m6Anet data.site_proba {path:?}:{line}: expected {} fields, found {}",
                SITE_PROBABILITY_COLUMNS.len(),
                record.len()
            );
        }
        let transcript_id = required(&record, 0, "transcript_id")?;
        validate_token("transcript_id", transcript_id)
            .with_context(|| format!("validate m6Anet data.site_proba {path:?}:{line}"))?;
        let pos0 = parse_u64(&record, 1, "transcript_position")
            .with_context(|| format!("validate m6Anet data.site_proba {path:?}:{line}"))?;
        let n_reads = parse_u64(&record, 2, "n_reads")
            .with_context(|| format!("validate m6Anet data.site_proba {path:?}:{line}"))?;
        let site_probability = parse_unit_interval(&record, 3, "probability_modified")
            .with_context(|| format!("validate m6Anet data.site_proba {path:?}:{line}"))?;
        let kmer = required(&record, 4, "kmer")?;
        validate_token("kmer", kmer)
            .with_context(|| format!("validate m6Anet data.site_proba {path:?}:{line}"))?;
        if !crate::modification::types::sequence_matches_iupac_motif(kmer, candidate_rule) {
            anyhow::bail!(
                "m6Anet data.site_proba {path:?}:{line}: kmer {kmer:?} does not match candidate_rule {candidate_rule:?}"
            );
        }
        let mod_ratio = parse_unit_interval(&record, 5, "mod_ratio")
            .with_context(|| format!("validate m6Anet data.site_proba {path:?}:{line}"))?;

        let site = TranscriptSite {
            transcript_id: transcript_id.to_owned(),
            pos0,
        };
        if !seen.insert(site.clone()) {
            anyhow::bail!(
                "m6Anet data.site_proba {path:?}:{line}: duplicate site {:?}:{}",
                site.transcript_id,
                site.pos0
            );
        }
        let probabilities = observed_site_probabilities.get(&site).with_context(|| {
            format!(
                "m6Anet data.site_proba {path:?}:{line}: site {:?}:{} is absent from data.indiv_proba",
                site.transcript_id, site.pos0
            )
        })?;
        let observed_n_reads =
            u64::try_from(probabilities.len()).context("m6Anet site read count exceeds u64")?;
        if n_reads != observed_n_reads {
            anyhow::bail!(
                "m6Anet data.site_proba {path:?}:{line}: site {:?}:{} has n_reads={n_reads}, but data.indiv_proba contains {observed_n_reads} rows",
                site.transcript_id,
                site.pos0
            );
        }
        if let Some(threshold) = read_probability_threshold {
            let modified = probabilities
                .iter()
                .filter(|probability| **probability >= threshold)
                .count();
            let expected_ratio = modified as f64 / probabilities.len() as f64;
            if (expected_ratio - mod_ratio).abs() > 1e-12 {
                anyhow::bail!(
                    "m6Anet data.site_proba {path:?}:{line}: site {:?}:{} has mod_ratio={mod_ratio}, expected {expected_ratio} from data.indiv_proba at read threshold {threshold}",
                    site.transcript_id,
                    site.pos0
                );
            }
        }
        total_reads = total_reads
            .checked_add(n_reads)
            .context("m6Anet data.site_proba total n_reads overflow")?;
        if site_probability >= SITE_PROBABILITY_THRESHOLD {
            sites_at_or_above_probability_threshold += 1;
        }
    }
    if seen.is_empty() {
        anyhow::bail!("m6Anet data.site_proba {path:?} contains no sites");
    }
    for site in observed_site_probabilities.keys() {
        if !seen.contains(site) {
            anyhow::bail!(
                "m6Anet data.site_proba {path:?} is missing retained site {:?}:{} from data.indiv_proba",
                site.transcript_id,
                site.pos0
            );
        }
    }

    Ok(M6anetSiteProbabilityQc {
        sites: seen.len(),
        total_reads,
        sites_at_or_above_probability_threshold,
        site_probability_threshold: SITE_PROBABILITY_THRESHOLD,
        read_probability_threshold,
    })
}

fn source_path(path: &Path) -> String {
    path.to_string_lossy().into_owned()
}

/// Import m6Anet RNA002 read probabilities into normalized genomic observations.
///
/// `read_map_path` is a strict tab-delimited file with exactly the header
/// `read_index<TAB>read_id`. Both the probability input and annotation may be
/// plain text or gzip-compressed. Transcript IDs are matched exactly, including
/// version suffixes, and `transcript_position` is interpreted as the zero-based
/// target base in 5'-to-3' spliced transcript space.
pub fn import_m6anet(
    data_indiv_proba_path: &Path,
    read_map_path: &Path,
    annotation_path: &Path,
    data_info_path: Option<&Path>,
    options: &M6anetImportOptions,
) -> anyhow::Result<M6anetImportResult> {
    import_m6anet_with_site_probability(
        data_indiv_proba_path,
        read_map_path,
        annotation_path,
        data_info_path,
        None,
        options,
    )
}

/// Import m6Anet observations and optionally cross-check `data.site_proba.csv`.
///
/// The site table contributes QC only. Read-to-isoform observations always come
/// exclusively from `data.indiv_proba.csv` and the explicit read-index map.
pub fn import_m6anet_with_site_probability(
    data_indiv_proba_path: &Path,
    read_map_path: &Path,
    annotation_path: &Path,
    data_info_path: Option<&Path>,
    data_site_probability_path: Option<&Path>,
    options: &M6anetImportOptions,
) -> anyhow::Result<M6anetImportResult> {
    let candidate_rule = validate_options(options).context("validate m6Anet import options")?;
    let read_map = read_read_map(read_map_path, &options.sample)?;
    let mut source_transcript_ids = read_source_transcript_ids(
        data_indiv_proba_path,
        &INDIVIDUAL_COLUMNS,
        "m6Anet data.indiv_proba",
    )?;
    if let Some(path) = data_info_path {
        source_transcript_ids.extend(read_source_transcript_ids(
            path,
            &DATA_INFO_COLUMNS,
            "m6Anet data.info",
        )?);
    }
    let source_transcript_bases = source_transcript_ids
        .iter()
        .map(|id| unversioned_transcript_id(id).to_owned())
        .collect::<BTreeSet<_>>();
    let transcripts = read_annotation_transcripts_where(
        annotation_path,
        &options.annotation_options,
        false,
        |annotation_id| {
            source_transcript_ids.contains(annotation_id)
                || source_transcript_bases.contains(unversioned_transcript_id(annotation_id))
        },
    )
    .with_context(|| format!("read m6Anet projection annotation {annotation_path:?}"))?;
    let transcript_index = TranscriptIndex::new(&transcripts)
        .context("index m6Anet projection annotation transcripts")?;

    let input = open_maybe_gzip(data_indiv_proba_path)
        .with_context(|| format!("open m6Anet data.indiv_proba {data_indiv_proba_path:?}"))?;
    let mut reader = ReaderBuilder::new()
        .delimiter(b',')
        .has_headers(true)
        .flexible(false)
        .from_reader(input);
    exact_header(
        &mut reader,
        &INDIVIDUAL_COLUMNS,
        "m6Anet data.indiv_proba",
        data_indiv_proba_path,
    )?;

    let mut input_rows = 0usize;
    let mut duplicate_exact = 0usize;
    let mut observations: BTreeMap<ModObservationKey, (usize, ModObservation)> = BTreeMap::new();
    let mut used_read_indices = BTreeSet::new();
    let mut used_transcripts = BTreeSet::new();
    let mut observed_site_reads: BTreeMap<TranscriptSite, u64> = BTreeMap::new();
    let mut observed_site_probabilities: BTreeMap<TranscriptSite, Vec<f64>> = BTreeMap::new();

    for (record_index, result) in reader.records().enumerate() {
        let line = record_index + 2;
        let record = result.with_context(|| {
            format!("parse m6Anet data.indiv_proba {data_indiv_proba_path:?}:{line}")
        })?;
        if record.len() != INDIVIDUAL_COLUMNS.len() {
            anyhow::bail!(
                "m6Anet data.indiv_proba {data_indiv_proba_path:?}:{line}: expected {} fields, found {}",
                INDIVIDUAL_COLUMNS.len(),
                record.len()
            );
        }
        input_rows += 1;

        let transcript_id = required(&record, 0, "transcript_id")?;
        validate_token("transcript_id", transcript_id).with_context(|| {
            format!("validate m6Anet data.indiv_proba {data_indiv_proba_path:?}:{line}")
        })?;
        let source_pos0 = parse_u64(&record, 1, "transcript_position").with_context(|| {
            format!("validate m6Anet data.indiv_proba {data_indiv_proba_path:?}:{line}")
        })?;
        let read_index = required(&record, 2, "read_index")?;
        validate_token("read_index", read_index).with_context(|| {
            format!("validate m6Anet data.indiv_proba {data_indiv_proba_path:?}:{line}")
        })?;
        let probability = parse_probability(&record, 3).with_context(|| {
            format!("validate m6Anet data.indiv_proba {data_indiv_proba_path:?}:{line}")
        })?;

        let mapped_read = read_map.get(read_index).with_context(|| {
            format!(
                "m6Anet data.indiv_proba {data_indiv_proba_path:?}:{line}: read_index {read_index} is missing from explicit read map {read_map_path:?}"
            )
        })?;
        let transcript = transcript_index
            .resolve(&transcripts, transcript_id)
            .with_context(|| {
                format!("validate m6Anet data.indiv_proba {data_indiv_proba_path:?}:{line}")
            })?;
        let genomic_pos0 = project_transcript_position(transcript, transcript_id, source_pos0)
            .with_context(|| {
                format!("validate m6Anet data.indiv_proba {data_indiv_proba_path:?}:{line}")
            })?;

        let observation = ModObservation {
            key: ModObservationKey {
                assay_id: options.assay_id.clone(),
                sample: options.sample.clone(),
                read_id: mapped_read.read_id.clone(),
                site: ModSiteKey {
                    chrom: transcript.chrom.clone(),
                    pos0: genomic_pos0,
                    strand: transcript.strand,
                    mod_code: M6A_MOD_CODE.to_owned(),
                },
            },
            probability: Some(probability),
            observation_state: ObservationState::ExplicitProbability,
            context: None,
            source_transcript_id: Some(transcript_id.to_owned()),
            source_pos0: Some(source_pos0),
        };
        observation
            .validate()
            .map_err(anyhow::Error::msg)
            .with_context(|| {
                format!(
                    "normalize m6Anet data.indiv_proba {data_indiv_proba_path:?}:{line} using read map line {}",
                    mapped_read.line
                )
            })?;

        used_read_indices.insert(read_index.to_owned());
        used_transcripts.insert(transcript_id.to_owned());
        let site = TranscriptSite {
            transcript_id: transcript_id.to_owned(),
            pos0: source_pos0,
        };

        match observations.get(&observation.key) {
            Some((_, existing)) if existing == &observation => duplicate_exact += 1,
            Some((first_line, _)) => {
                anyhow::bail!(
                    "m6Anet data.indiv_proba {data_indiv_proba_path:?}:{line} conflicts with duplicate normalized key first seen at line {first_line}: {:?}",
                    observation.key
                );
            }
            None => {
                let count = observed_site_reads.entry(site.clone()).or_default();
                *count = count
                    .checked_add(1)
                    .context("m6Anet retained site read count overflow")?;
                observed_site_probabilities
                    .entry(site)
                    .or_default()
                    .push(probability);
                observations.insert(observation.key.clone(), (line, observation));
            }
        }
    }
    if input_rows == 0 {
        anyhow::bail!("m6Anet data.indiv_proba {data_indiv_proba_path:?} contains no observations");
    }

    let data_info = data_info_path
        .map(|path| {
            audit_data_info(
                path,
                &observed_site_reads,
                options.minimum_reads,
                &transcripts,
                &transcript_index,
            )
        })
        .transpose()?;
    let site_probability = data_site_probability_path
        .map(|path| {
            audit_site_probability(
                path,
                &observed_site_probabilities,
                options.read_probability_threshold,
                &candidate_rule,
            )
        })
        .transpose()?;

    let source_site_filter = match (data_info.is_some(), options.minimum_reads) {
        (true, Some(minimum_reads)) => {
            format!("m6anet_data_info_n_reads>={minimum_reads}")
        }
        (true, None) => "m6anet_data_info_retained_sites".to_owned(),
        (false, Some(minimum_reads)) => {
            format!("m6anet_min_reads={minimum_reads}_unverified_without_data_info")
        }
        (false, None) => "m6anet_retained_sites_unverified".to_owned(),
    };
    let mut source_files = vec![
        source_path(data_indiv_proba_path),
        source_path(read_map_path),
        source_path(annotation_path),
    ];
    if let Some(path) = data_info_path {
        source_files.push(source_path(path));
    }
    if let Some(path) = data_site_probability_path {
        source_files.push(source_path(path));
    }
    let metadata = AssayMetadata {
        schema_version: MODIFICATION_SCHEMA_VERSION,
        assay_id: options.assay_id.clone(),
        caller: "m6anet".to_owned(),
        caller_version: options.caller_version.clone(),
        model_id: options.model_id.clone(),
        chemistry: "RNA002".to_owned(),
        candidate_rule,
        source_emission_threshold: None,
        source_site_filter,
        candidate_observations_complete: data_info.is_some(),
        implicit_skip_policy: ImplicitSkipPolicy::NotApplicable,
        coordinate_source:
            "m6anet_transcript_position_0_based_target_base_projected_via_annotation".to_owned(),
        read_id_mapping: "explicit_read_index_read_id_map".to_owned(),
        source_files,
    };
    metadata
        .validate()
        .map_err(anyhow::Error::msg)
        .context("validate imported m6Anet assay metadata")?;

    let observations: Vec<ModObservation> = observations
        .into_values()
        .map(|(_, observation)| observation)
        .collect();
    let qc = M6anetImportQc {
        input_rows,
        unique_observations: observations.len(),
        duplicate_exact,
        read_map_entries: read_map.len(),
        read_map_entries_used: used_read_indices.len(),
        source_transcripts: used_transcripts.len(),
        projection_transcripts_loaded: transcripts.len(),
        source_sites: observed_site_reads.len(),
        data_info,
        site_probability,
    };
    Ok(M6anetImportResult {
        observations,
        metadata,
        qc,
    })
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::io::Write;
    use std::path::{Path, PathBuf};
    use std::time::{SystemTime, UNIX_EPOCH};

    use flate2::write::GzEncoder;
    use flate2::Compression;

    use super::*;

    fn temp_dir(label: &str) -> PathBuf {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("time")
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "trackcluster-m6anet-{label}-{}-{nonce}",
            std::process::id()
        ));
        fs::create_dir(&path).unwrap();
        path
    }

    fn write_gzip(path: &Path, contents: &str) {
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(contents.as_bytes()).unwrap();
        fs::write(path, encoder.finish().unwrap()).unwrap();
    }

    fn annotation() -> &'static str {
        concat!(
            "chr1\ttest\texon\t101\t103\t.\t+\t.\tgene_id \"G1\"; transcript_id \"TXPLUS.1\";\n",
            "chr1\ttest\texon\t201\t203\t.\t+\t.\tgene_id \"G1\"; transcript_id \"TXPLUS.1\";\n",
            "chr2\ttest\texon\t1001\t1003\t.\t-\t.\tgene_id \"G2\"; transcript_id \"TXMINUS.2\";\n",
            "chr2\ttest\texon\t2001\t2003\t.\t-\t.\tgene_id \"G2\"; transcript_id \"TXMINUS.2\";\n",
        )
    }

    fn options() -> M6anetImportOptions {
        let mut options = M6anetImportOptions::new("S1", "m6anet_rna002", "HCT116_RNA002");
        options.caller_version = "2.1.0".to_owned();
        options
    }

    fn write_common(root: &Path, annotation_text: &str, read_map_text: &str) -> (PathBuf, PathBuf) {
        let annotation_path = root.join("annotation.gtf");
        let read_map_path = root.join("read-map.tsv");
        fs::write(&annotation_path, annotation_text).unwrap();
        fs::write(&read_map_path, read_map_text).unwrap();
        (annotation_path, read_map_path)
    }

    fn error_text(result: anyhow::Result<M6anetImportResult>) -> String {
        format!("{:#}", result.unwrap_err())
    }

    #[test]
    fn candidate_rule_is_normalized_and_requires_an_a_centered_iupac_motif() {
        let mut options = options();
        options.candidate_rule = "drach".to_owned();
        assert_eq!(validate_options(&options).unwrap(), "DRACH");

        for invalid in ["foo", "CCCCC", "DRA?H", "NNNN"] {
            options.candidate_rule = invalid.to_owned();
            assert!(
                validate_options(&options).is_err(),
                "{invalid:?} unexpectedly passed"
            );
        }
    }

    #[test]
    fn site_probability_kmer_must_match_the_candidate_rule() {
        let root = temp_dir("site-kmer");
        let site_probability_path = root.join("data.site_proba.csv");
        fs::write(
            &site_probability_path,
            concat!(
                "transcript_id,transcript_position,n_reads,probability_modified,kmer,mod_ratio\n",
                "TXPLUS.1,3,1,0.9,CCCCC,1\n",
            ),
        )
        .unwrap();
        let observed = BTreeMap::from([(
            TranscriptSite {
                transcript_id: "TXPLUS.1".to_owned(),
                pos0: 3,
            },
            vec![0.9],
        )]);

        let error = audit_site_probability(&site_probability_path, &observed, Some(0.5), "DRACH")
            .unwrap_err();
        assert!(format!("{error:#}").contains("does not match candidate_rule"));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn imports_plain_and_gzip_with_strand_aware_projection_and_prefix_audit() {
        let root = temp_dir("projection");
        let annotation_path = root.join("annotation.gtf.gz");
        write_gzip(&annotation_path, annotation());
        let read_map_path = root.join("read-map.tsv");
        fs::write(
            &read_map_path,
            concat!(
                "read_index\tread_id\n",
                "1\tread-plus\n",
                "2\tS1::read-minus\n",
                "3\tunused\n",
            ),
        )
        .unwrap();
        let individual = concat!(
            "transcript_id,transcript_position,read_index,probability_modified\n",
            "TXPLUS.1,3,1,0.25\n",
            "TXMINUS.2,0,2,0.75\n",
        );
        let plain_path = root.join("data.indiv_proba.csv");
        let gzip_path = root.join("data.indiv_proba.csv.gz");
        fs::write(&plain_path, individual).unwrap();
        write_gzip(&gzip_path, individual);
        let data_info_path = root.join("data.info");
        fs::write(
            &data_info_path,
            concat!(
                "transcript_id,transcript_position,start,end,n_reads\n",
                "TXPLUS.1,3,0,10,1\n",
                "TXMINUS.2,0,10,20,1\n",
            ),
        )
        .unwrap();
        let mut options = options();
        options.minimum_reads = Some(1);

        let plain = import_m6anet(
            &plain_path,
            &read_map_path,
            &annotation_path,
            Some(&data_info_path),
            &options,
        )
        .unwrap();
        let gzip = import_m6anet(
            &gzip_path,
            &read_map_path,
            &annotation_path,
            Some(&data_info_path),
            &options,
        )
        .unwrap();
        assert_eq!(plain.observations, gzip.observations);
        assert_eq!(plain.observations.len(), 2);

        let plus = plain
            .observations
            .iter()
            .find(|observation| observation.source_transcript_id.as_deref() == Some("TXPLUS.1"))
            .unwrap();
        assert_eq!(plus.key.read_id, "S1::read-plus");
        assert_eq!(plus.key.site.chrom, "chr1");
        assert_eq!(plus.key.site.pos0, 200);
        assert_eq!(plus.key.site.strand, Strand::Plus);
        assert_eq!(plus.key.site.mod_code, "A+a");
        assert_eq!(plus.source_pos0, Some(3));

        let minus = plain
            .observations
            .iter()
            .find(|observation| observation.source_transcript_id.as_deref() == Some("TXMINUS.2"))
            .unwrap();
        assert_eq!(minus.key.read_id, "S1::read-minus");
        assert_eq!(minus.key.site.chrom, "chr2");
        assert_eq!(minus.key.site.pos0, 2002);
        assert_eq!(minus.key.site.strand, Strand::Minus);
        assert_eq!(minus.source_pos0, Some(0));

        assert_eq!(plain.qc.input_rows, 2);
        assert_eq!(plain.qc.read_map_entries, 3);
        assert_eq!(plain.qc.read_map_entries_used, 2);
        assert_eq!(plain.qc.source_transcripts, 2);
        assert_eq!(plain.qc.projection_transcripts_loaded, 2);
        assert_eq!(plain.qc.source_sites, 2);
        assert_eq!(plain.qc.data_info.as_ref().unwrap().sites, 2);
        assert!(plain.metadata.candidate_observations_complete);
        assert_eq!(plain.metadata.chemistry, "RNA002");
        assert_eq!(
            plain.metadata.coordinate_source,
            "m6anet_transcript_position_0_based_target_base_projected_via_annotation"
        );
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn read_map_requires_exact_header_complete_mapping_and_matching_prefix() {
        let root = temp_dir("read-map-errors");
        let individual_path = root.join("data.indiv_proba.csv");
        fs::write(
            &individual_path,
            concat!(
                "transcript_id,transcript_position,read_index,probability_modified\n",
                "TXPLUS.1,0,1,0.5\n",
            ),
        )
        .unwrap();
        let (annotation_path, read_map_path) =
            write_common(&root, annotation(), "read_index\tname\n1\tread1\n");
        let error = error_text(import_m6anet(
            &individual_path,
            &read_map_path,
            &annotation_path,
            None,
            &options(),
        ));
        assert!(error.contains("header mismatch"));

        fs::write(&read_map_path, "read_index\tread_id\n2\tread2\n").unwrap();
        let error = error_text(import_m6anet(
            &individual_path,
            &read_map_path,
            &annotation_path,
            None,
            &options(),
        ));
        assert!(error.contains("read_index 1 is missing"));

        fs::write(&read_map_path, "read_index\tread_id\n1\tS2::read1\n").unwrap();
        let error = error_text(import_m6anet(
            &individual_path,
            &read_map_path,
            &annotation_path,
            None,
            &options(),
        ));
        assert!(error.contains("sample prefix \"S2\", expected \"S1\""));

        fs::write(&read_map_path, "read_index\tread_id\n").unwrap();
        let error = error_text(import_m6anet(
            &individual_path,
            &read_map_path,
            &annotation_path,
            None,
            &options(),
        ));
        assert!(error.contains("contains no mappings"));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn read_index_is_an_opaque_token_and_accepts_replicate_suffixes() {
        let root = temp_dir("opaque-read-index");
        let (annotation_path, read_map_path) = write_common(
            &root,
            annotation(),
            "read_index\tread_id\n966210_0\tread1\n",
        );
        let individual_path = root.join("data.indiv_proba.csv");
        fs::write(
            &individual_path,
            concat!(
                "transcript_id,transcript_position,read_index,probability_modified\n",
                "TXPLUS.1,0,966210_0,0.5\n",
            ),
        )
        .unwrap();

        let imported = import_m6anet(
            &individual_path,
            &read_map_path,
            &annotation_path,
            None,
            &options(),
        )
        .unwrap();
        assert_eq!(imported.observations.len(), 1);
        assert_eq!(imported.observations[0].key.read_id, "S1::read1");
        assert_eq!(imported.qc.read_map_entries_used, 1);
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn transcript_lookup_requires_exact_id_and_version() {
        let root = temp_dir("transcript-identity");
        let (annotation_path, read_map_path) =
            write_common(&root, annotation(), "read_index\tread_id\n1\tread1\n");
        let individual_path = root.join("data.indiv_proba.csv");
        fs::write(
            &individual_path,
            concat!(
                "transcript_id,transcript_position,read_index,probability_modified\n",
                "TXPLUS.2,0,1,0.5\n",
            ),
        )
        .unwrap();
        let error = error_text(import_m6anet(
            &individual_path,
            &read_map_path,
            &annotation_path,
            None,
            &options(),
        ));
        assert!(error.contains("transcript version mismatch"));
        assert!(error.contains("TXPLUS.1"));

        fs::write(
            &individual_path,
            concat!(
                "transcript_id,transcript_position,read_index,probability_modified\n",
                "MISSING.1,0,1,0.5\n",
            ),
        )
        .unwrap();
        let error = error_text(import_m6anet(
            &individual_path,
            &read_map_path,
            &annotation_path,
            None,
            &options(),
        ));
        assert!(error.contains("is absent from the annotation"));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn projection_rejects_out_of_bounds_and_unknown_strand() {
        let root = temp_dir("projection-errors");
        let (annotation_path, read_map_path) =
            write_common(&root, annotation(), "read_index\tread_id\n1\tread1\n");
        let individual_path = root.join("data.indiv_proba.csv");
        fs::write(
            &individual_path,
            concat!(
                "transcript_id,transcript_position,read_index,probability_modified\n",
                "TXPLUS.1,6,1,0.5\n",
            ),
        )
        .unwrap();
        let error = error_text(import_m6anet(
            &individual_path,
            &read_map_path,
            &annotation_path,
            None,
            &options(),
        ));
        assert!(error.contains("out of bounds"));
        assert!(error.contains("spliced length 6"));

        fs::write(
            &annotation_path,
            "chr1\ttest\texon\t101\t103\t.\t.\t.\tgene_id \"G1\"; transcript_id \"UNKNOWN.1\";\n",
        )
        .unwrap();
        fs::write(
            &individual_path,
            concat!(
                "transcript_id,transcript_position,read_index,probability_modified\n",
                "UNKNOWN.1,0,1,0.5\n",
            ),
        )
        .unwrap();
        let error = error_text(import_m6anet(
            &individual_path,
            &read_map_path,
            &annotation_path,
            None,
            &options(),
        ));
        assert!(error.contains("has unknown strand"));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn probability_and_coordinate_fields_are_strict() {
        let root = temp_dir("numeric-errors");
        let (annotation_path, read_map_path) =
            write_common(&root, annotation(), "read_index\tread_id\n1\tread1\n");
        let individual_path = root.join("data.indiv_proba.csv");
        for value in ["NaN", "inf", "-0.1", "1.1", "not-a-number"] {
            fs::write(
                &individual_path,
                format!(
                    "transcript_id,transcript_position,read_index,probability_modified\nTXPLUS.1,0,1,{value}\n"
                ),
            )
            .unwrap();
            let error = error_text(import_m6anet(
                &individual_path,
                &read_map_path,
                &annotation_path,
                None,
                &options(),
            ));
            assert!(error.contains("probability_modified"), "{error}");
        }

        fs::write(
            &individual_path,
            concat!(
                "transcript_id,transcript_position,read_index,probability_modified\n",
                "TXPLUS.1,not-an-offset,1,0.5\n",
            ),
        )
        .unwrap();
        let error = error_text(import_m6anet(
            &individual_path,
            &read_map_path,
            &annotation_path,
            None,
            &options(),
        ));
        assert!(error.contains("invalid unsigned integer for transcript_position"));

        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn exact_duplicates_fold_but_conflicting_duplicates_fail() {
        let root = temp_dir("duplicates");
        let (annotation_path, read_map_path) =
            write_common(&root, annotation(), "read_index\tread_id\n1\tread1\n");
        let individual_path = root.join("data.indiv_proba.csv");
        fs::write(
            &individual_path,
            concat!(
                "transcript_id,transcript_position,read_index,probability_modified\n",
                "TXPLUS.1,0,1,0.1\n",
                "TXPLUS.1,0,1,0.1\n",
            ),
        )
        .unwrap();
        let imported = import_m6anet(
            &individual_path,
            &read_map_path,
            &annotation_path,
            None,
            &options(),
        )
        .unwrap();
        assert_eq!(imported.qc.input_rows, 2);
        assert_eq!(imported.qc.duplicate_exact, 1);
        assert_eq!(imported.qc.unique_observations, 1);
        assert!(!imported.metadata.candidate_observations_complete);

        let data_info_path = root.join("data.info");
        fs::write(
            &data_info_path,
            concat!(
                "transcript_id,transcript_position,start,end,n_reads\n",
                "TXPLUS.1,0,0,10,1\n",
            ),
        )
        .unwrap();
        let audited = import_m6anet(
            &individual_path,
            &read_map_path,
            &annotation_path,
            Some(&data_info_path),
            &options(),
        )
        .unwrap();
        assert_eq!(audited.qc.data_info.unwrap().retained_reads, 1);

        fs::write(
            &individual_path,
            concat!(
                "transcript_id,transcript_position,read_index,probability_modified\n",
                "TXPLUS.1,0,1,0.1\n",
                "TXPLUS.1,0,1,0.2\n",
            ),
        )
        .unwrap();
        let error = error_text(import_m6anet(
            &individual_path,
            &read_map_path,
            &annotation_path,
            None,
            &options(),
        ));
        assert!(error.contains("conflicts with duplicate normalized key"));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn data_info_audits_n_reads_and_exact_retention_threshold() {
        let root = temp_dir("data-info");
        let (annotation_path, read_map_path) = write_common(
            &root,
            annotation(),
            concat!("read_index\tread_id\n", "1\tread1\n", "2\tread2\n",),
        );
        let individual_path = root.join("data.indiv_proba.csv");
        fs::write(
            &individual_path,
            concat!(
                "transcript_id,transcript_position,read_index,probability_modified\n",
                "TXPLUS.1,3,1,0.1\n",
                "TXPLUS.1,3,2,0.2\n",
            ),
        )
        .unwrap();
        let data_info_path = root.join("data.info");
        fs::write(
            &data_info_path,
            concat!(
                "transcript_id,transcript_position,start,end,n_reads\n",
                "TXPLUS.1,3,0,20,2\n",
                "TXPLUS.1,4,20,30,1\n",
            ),
        )
        .unwrap();
        let mut options = options();
        options.minimum_reads = Some(2);
        let imported = import_m6anet(
            &individual_path,
            &read_map_path,
            &annotation_path,
            Some(&data_info_path),
            &options,
        )
        .unwrap();
        assert_eq!(
            imported.qc.data_info,
            Some(M6anetDataInfoQc {
                sites: 2,
                retained_sites: 1,
                filtered_sites: 1,
                total_reads: 3,
                retained_reads: 2,
                minimum_reads: Some(2),
            })
        );
        assert_eq!(
            imported.metadata.source_site_filter,
            "m6anet_data_info_n_reads>=2"
        );

        fs::write(
            &data_info_path,
            concat!(
                "transcript_id,transcript_position,start,end,n_reads\n",
                "TXPLUS.1,3,0,20,3\n",
                "TXPLUS.1,4,20,30,1\n",
            ),
        )
        .unwrap();
        let error = error_text(import_m6anet(
            &individual_path,
            &read_map_path,
            &annotation_path,
            Some(&data_info_path),
            &options,
        ));
        assert!(error.contains("data.indiv_proba contains 2 rows"));

        fs::write(
            &individual_path,
            concat!(
                "transcript_id,transcript_position,read_index,probability_modified\n",
                "TXPLUS.1,3,1,0.1\n",
            ),
        )
        .unwrap();
        fs::write(
            &data_info_path,
            concat!(
                "transcript_id,transcript_position,start,end,n_reads\n",
                "TXPLUS.1,3,0,20,1\n",
            ),
        )
        .unwrap();
        let error = error_text(import_m6anet(
            &individual_path,
            &read_map_path,
            &annotation_path,
            Some(&data_info_path),
            &options,
        ));
        assert!(error.contains("retained-site mismatch"));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn known_rna002_models_have_distinct_official_read_threshold_presets() {
        assert_eq!(
            read_probability_threshold_preset("HCT116_RNA002"),
            Some(HCT116_RNA002_READ_THRESHOLD)
        );
        assert_eq!(
            read_probability_threshold_preset("arabidopsis_RNA002"),
            Some(ARABIDOPSIS_RNA002_READ_THRESHOLD)
        );
        assert_ne!(
            HCT116_RNA002_READ_THRESHOLD,
            ARABIDOPSIS_RNA002_READ_THRESHOLD
        );
        assert_eq!(read_probability_threshold_preset("custom_model"), None);
    }
}
