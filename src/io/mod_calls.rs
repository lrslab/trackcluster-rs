//! Strict normalized modification observation and assay metadata I/O.

use std::collections::BTreeMap;
use std::fs::File;
use std::io::Write;
use std::path::Path;

use anyhow::Context;
use csv::{ReaderBuilder, StringRecord, WriterBuilder};

use crate::model::Strand;
use crate::modification::{
    AssayMetadata, ModObservation, ModObservationKey, ModSiteKey, ObservationState,
};

/// V1 normalized observation TSV columns in their required order.
pub const OBSERVATION_COLUMNS: [&str; 12] = [
    "assay_id",
    "sample",
    "read_id",
    "chrom",
    "pos0",
    "strand",
    "mod_code",
    "probability",
    "observation_state",
    "context",
    "source_transcript_id",
    "source_pos0",
];

/// Parsed normalized observations and deterministic duplicate audit counts.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ObservationReadResult {
    /// Unique validated observations in deterministic key order.
    pub observations: Vec<ModObservation>,
    /// Physical data rows parsed, including exact duplicates.
    pub input_rows: usize,
    /// Exact duplicate rows removed.
    pub duplicate_exact: usize,
}

/// Audit counts from a canonical streaming observation scan.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ObservationScanResult {
    /// Physical data rows parsed, including exact duplicates.
    pub input_rows: usize,
    /// Unique validated observations passed to the visitor.
    pub valid_rows: usize,
    /// Adjacent exact duplicate rows omitted.
    pub duplicate_exact: usize,
}

fn required<'a>(record: &'a StringRecord, index: usize, name: &str) -> anyhow::Result<&'a str> {
    let value = record.get(index).unwrap_or_default();
    if value.is_empty() {
        anyhow::bail!("{name} must not be empty");
    }
    Ok(value)
}

fn parse_optional(value: &str, name: &str) -> anyhow::Result<Option<String>> {
    match value {
        "NA" => Ok(None),
        "" => anyhow::bail!("{name} must use NA for a missing value"),
        value => Ok(Some(value.to_owned())),
    }
}

fn parse_observation(record: &StringRecord) -> anyhow::Result<ModObservation> {
    if record.len() != OBSERVATION_COLUMNS.len() {
        anyhow::bail!(
            "expected {} fields, found {}",
            OBSERVATION_COLUMNS.len(),
            record.len()
        );
    }

    let strand = Strand::try_from(required(record, 5, "strand")?)?;
    if strand == Strand::Unknown {
        anyhow::bail!("normalized modification strand must be '+' or '-'");
    }
    let probability_field = required(record, 7, "probability")?;
    let probability = if probability_field == "NA" {
        None
    } else {
        Some(
            probability_field
                .parse::<f64>()
                .with_context(|| format!("invalid probability {probability_field:?}"))?,
        )
    };
    let source_pos_field = required(record, 11, "source_pos0")?;
    let source_pos0 = if source_pos_field == "NA" {
        None
    } else {
        Some(
            source_pos_field
                .parse::<u64>()
                .with_context(|| format!("invalid source_pos0 {source_pos_field:?}"))?,
        )
    };

    let observation = ModObservation {
        key: ModObservationKey {
            assay_id: required(record, 0, "assay_id")?.to_owned(),
            sample: required(record, 1, "sample")?.to_owned(),
            read_id: required(record, 2, "read_id")?.to_owned(),
            site: ModSiteKey {
                chrom: required(record, 3, "chrom")?.to_owned(),
                pos0: required(record, 4, "pos0")?
                    .parse::<u32>()
                    .with_context(|| format!("invalid pos0 {:?}", record.get(4)))?,
                strand,
                mod_code: required(record, 6, "mod_code")?.to_owned(),
            },
        },
        probability,
        observation_state: required(record, 8, "observation_state")?
            .parse::<ObservationState>()
            .map_err(anyhow::Error::msg)?,
        context: parse_optional(required(record, 9, "context")?, "context")?,
        source_transcript_id: parse_optional(
            required(record, 10, "source_transcript_id")?,
            "source_transcript_id",
        )?,
        source_pos0,
    };
    observation.validate().map_err(anyhow::Error::msg)?;
    Ok(observation)
}

/// Read, validate, sort, and exactly deduplicate a normalized observation TSV.
pub fn read_observations_tsv(path: &Path) -> anyhow::Result<ObservationReadResult> {
    let file = File::open(path).with_context(|| format!("open observations {path:?}"))?;
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .flexible(false)
        .from_reader(file);

    let expected = StringRecord::from(OBSERVATION_COLUMNS.to_vec());
    let actual = reader
        .headers()
        .with_context(|| format!("read observations header {path:?}"))?;
    if actual != &expected {
        anyhow::bail!(
            "observations {path:?} header mismatch: expected {:?}, found {:?}",
            OBSERVATION_COLUMNS,
            actual.iter().collect::<Vec<_>>()
        );
    }

    let mut input_rows = 0usize;
    let mut duplicate_exact = 0usize;
    let mut by_key: BTreeMap<ModObservationKey, (usize, ModObservation)> = BTreeMap::new();
    for (record_index, result) in reader.records().enumerate() {
        let line = record_index + 2;
        let record = result.with_context(|| format!("parse observations {path:?}:{line}"))?;
        input_rows += 1;
        let observation = parse_observation(&record)
            .with_context(|| format!("validate observations {path:?}:{line}"))?;
        match by_key.get(&observation.key) {
            Some((_, existing)) if existing == &observation => duplicate_exact += 1,
            Some((first_line, _)) => {
                anyhow::bail!(
                    "observations {path:?}:{line} conflicts with duplicate key first seen at line {first_line}: {:?}",
                    observation.key
                );
            }
            None => {
                by_key.insert(observation.key.clone(), (line, observation));
            }
        }
    }

    Ok(ObservationReadResult {
        observations: by_key.into_values().map(|(_, value)| value).collect(),
        input_rows,
        duplicate_exact,
    })
}

/// Stream a canonically sorted observation TSV without materializing it.
///
/// Importer outputs are sorted by the full observation key. This scanner
/// verifies that invariant, folds adjacent exact duplicates, and rejects
/// conflicting duplicates before invoking `visit`.
pub fn scan_canonical_observations_tsv<F>(
    path: &Path,
    mut visit: F,
) -> anyhow::Result<ObservationScanResult>
where
    F: FnMut(&ModObservation) -> anyhow::Result<()>,
{
    let file = File::open(path).with_context(|| format!("open observations {path:?}"))?;
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .flexible(false)
        .from_reader(file);

    let expected = StringRecord::from(OBSERVATION_COLUMNS.to_vec());
    let actual = reader
        .headers()
        .with_context(|| format!("read observations header {path:?}"))?;
    if actual != &expected {
        anyhow::bail!(
            "observations {path:?} header mismatch: expected {:?}, found {:?}",
            OBSERVATION_COLUMNS,
            actual.iter().collect::<Vec<_>>()
        );
    }

    let mut result = ObservationScanResult::default();
    let mut previous: Option<(usize, ModObservation)> = None;
    for (record_index, record) in reader.records().enumerate() {
        let line = record_index + 2;
        let record = record.with_context(|| format!("parse observations {path:?}:{line}"))?;
        result.input_rows += 1;
        let observation = parse_observation(&record)
            .with_context(|| format!("validate observations {path:?}:{line}"))?;
        if let Some((previous_line, previous_observation)) = previous.as_ref() {
            match observation.key.cmp(&previous_observation.key) {
                std::cmp::Ordering::Less => {
                    anyhow::bail!(
                        "observations {path:?}:{line} are not in canonical key order; \
                         previous key at line {previous_line} was {:?}, current key is {:?}",
                        previous_observation.key,
                        observation.key
                    );
                }
                std::cmp::Ordering::Equal if &observation == previous_observation => {
                    result.duplicate_exact += 1;
                    continue;
                }
                std::cmp::Ordering::Equal => {
                    anyhow::bail!(
                        "observations {path:?}:{line} conflict with duplicate key first seen at \
                         line {previous_line}: {:?}",
                        observation.key
                    );
                }
                std::cmp::Ordering::Greater => {}
            }
        }
        visit(&observation).with_context(|| format!("process observations {path:?}:{line}"))?;
        result.valid_rows += 1;
        previous = Some((line, observation));
    }
    Ok(result)
}

fn optional_field(value: Option<&str>) -> &str {
    value.unwrap_or("NA")
}

pub(crate) fn new_observation_writer<W: Write>(writer: W) -> anyhow::Result<csv::Writer<W>> {
    let mut output = WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_writer(writer);
    output.write_record(OBSERVATION_COLUMNS)?;
    Ok(output)
}

pub(crate) fn write_observation_record<W: Write>(
    output: &mut csv::Writer<W>,
    observation: &ModObservation,
) -> anyhow::Result<()> {
    observation.validate().map_err(anyhow::Error::msg)?;
    let probability = observation
        .probability
        .map(|value| value.to_string())
        .unwrap_or_else(|| "NA".to_owned());
    let pos0 = observation.key.site.pos0.to_string();
    let strand = observation.key.site.strand.as_char().to_string();
    let source_pos0 = observation
        .source_pos0
        .map(|value| value.to_string())
        .unwrap_or_else(|| "NA".to_owned());
    output.write_record([
        observation.key.assay_id.as_str(),
        observation.key.sample.as_str(),
        observation.key.read_id.as_str(),
        observation.key.site.chrom.as_str(),
        pos0.as_str(),
        strand.as_str(),
        observation.key.site.mod_code.as_str(),
        probability.as_str(),
        observation.observation_state.as_str(),
        optional_field(observation.context.as_deref()),
        optional_field(observation.source_transcript_id.as_deref()),
        source_pos0.as_str(),
    ])?;
    Ok(())
}

/// Write validated observations in deterministic key order.
pub fn write_observations_tsv_to_writer<W: Write>(
    writer: W,
    observations: &[ModObservation],
) -> anyhow::Result<()> {
    let mut sorted = observations.to_vec();
    sorted.sort_by(|left, right| left.key.cmp(&right.key));
    for observation in &sorted {
        observation.validate().map_err(anyhow::Error::msg)?;
    }
    for pair in sorted.windows(2) {
        if pair[0].key == pair[1].key {
            anyhow::bail!("cannot write duplicate observation key: {:?}", pair[0].key);
        }
    }

    let mut output = new_observation_writer(writer)?;
    for observation in &sorted {
        write_observation_record(&mut output, observation)?;
    }
    output.flush().context("flush normalized observations")?;
    Ok(())
}

/// Read and validate V1 assay provenance JSON.
pub fn read_assay_metadata(path: &Path) -> anyhow::Result<AssayMetadata> {
    let bytes = std::fs::read(path).with_context(|| format!("read assay metadata {path:?}"))?;
    let metadata: AssayMetadata =
        serde_json::from_slice(&bytes).with_context(|| format!("parse assay metadata {path:?}"))?;
    metadata
        .validate()
        .map_err(anyhow::Error::msg)
        .with_context(|| format!("validate assay metadata {path:?}"))?;
    Ok(metadata)
}

/// Write validated V1 assay provenance JSON with a terminating newline.
pub fn write_assay_metadata_to_writer<W: Write>(
    mut writer: W,
    metadata: &AssayMetadata,
) -> anyhow::Result<()> {
    metadata.validate().map_err(anyhow::Error::msg)?;
    serde_json::to_writer_pretty(&mut writer, metadata).context("serialize assay metadata")?;
    writer
        .write_all(b"\n")
        .context("terminate assay metadata JSON")?;
    writer.flush().context("flush assay metadata JSON")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::path::PathBuf;
    use std::time::{SystemTime, UNIX_EPOCH};

    use crate::modification::ImplicitSkipPolicy;

    use super::*;

    fn temp_file(label: &str) -> PathBuf {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("time")
            .as_nanos();
        std::env::temp_dir().join(format!(
            "trackcluster-mod-calls-{label}-{}-{nonce}",
            std::process::id()
        ))
    }

    fn row(probability: &str, state: &str) -> String {
        format!("rna004\tS1\tS1::read1\tchr1\t10\t+\tA+a\t{probability}\t{state}\tDRACH\tNA\tNA\n")
    }

    #[test]
    fn observation_tsv_round_trips_and_folds_exact_duplicates() {
        let path = temp_file("roundtrip.tsv");
        let input = format!(
            "{}\n{}{}",
            OBSERVATION_COLUMNS.join("\t"),
            row("0", "explicit_probability"),
            row("0", "explicit_probability")
        );
        fs::write(&path, input).unwrap();

        let parsed = read_observations_tsv(&path).unwrap();
        assert_eq!(parsed.input_rows, 2);
        assert_eq!(parsed.duplicate_exact, 1);
        assert_eq!(parsed.observations.len(), 1);

        let mut output = Vec::new();
        write_observations_tsv_to_writer(&mut output, &parsed.observations).unwrap();
        assert!(String::from_utf8(output)
            .unwrap()
            .contains("\t0\texplicit_probability\t"));
        let _ = fs::remove_file(path);
    }

    #[test]
    fn observation_tsv_rejects_conflicts_and_invalid_state_probability_pairs() {
        let conflict = temp_file("conflict.tsv");
        let input = format!(
            "{}\n{}{}",
            OBSERVATION_COLUMNS.join("\t"),
            row("0.1", "explicit_probability"),
            row("0.2", "explicit_probability")
        );
        fs::write(&conflict, input).unwrap();
        assert!(read_observations_tsv(&conflict)
            .unwrap_err()
            .to_string()
            .contains("conflicts with duplicate key"));

        let invalid = temp_file("invalid.tsv");
        fs::write(
            &invalid,
            format!(
                "{}\n{}",
                OBSERVATION_COLUMNS.join("\t"),
                row("0.2", "unknown")
            ),
        )
        .unwrap();
        assert!(read_observations_tsv(&invalid).is_err());
        let _ = fs::remove_file(conflict);
        let _ = fs::remove_file(invalid);
    }

    #[test]
    fn canonical_streaming_scan_folds_adjacent_duplicates_and_rejects_unsorted_keys() {
        let canonical = temp_file("canonical-scan.tsv");
        let read1 = row("0.1", "explicit_probability");
        let read2 = read1.replace("S1::read1", "S1::read2");
        fs::write(
            &canonical,
            format!(
                "{}\n{}{}{}",
                OBSERVATION_COLUMNS.join("\t"),
                read1,
                read1,
                read2
            ),
        )
        .unwrap();
        let mut visited = Vec::new();
        let result = scan_canonical_observations_tsv(&canonical, |observation| {
            visited.push(observation.key.read_id.clone());
            Ok(())
        })
        .unwrap();
        assert_eq!(result.input_rows, 3);
        assert_eq!(result.valid_rows, 2);
        assert_eq!(result.duplicate_exact, 1);
        assert_eq!(visited, ["S1::read1", "S1::read2"]);

        let unsorted = temp_file("unsorted-scan.tsv");
        fs::write(
            &unsorted,
            format!("{}\n{}{}", OBSERVATION_COLUMNS.join("\t"), read2, read1),
        )
        .unwrap();
        assert!(scan_canonical_observations_tsv(&unsorted, |_| Ok(()))
            .unwrap_err()
            .to_string()
            .contains("not in canonical key order"));
        let _ = fs::remove_file(canonical);
        let _ = fs::remove_file(unsorted);
    }

    #[test]
    fn assay_metadata_json_round_trips_and_rejects_unknown_fields() {
        let metadata = AssayMetadata {
            schema_version: 1,
            assay_id: "rna002".to_owned(),
            caller: "m6anet".to_owned(),
            caller_version: "2.1.0".to_owned(),
            model_id: "HCT116_RNA002".to_owned(),
            chemistry: "RNA002".to_owned(),
            candidate_rule: "DRACH".to_owned(),
            source_emission_threshold: None,
            source_site_filter: "min_reads=20".to_owned(),
            candidate_observations_complete: true,
            implicit_skip_policy: ImplicitSkipPolicy::NotApplicable,
            coordinate_source: "reference_transcript_projected_to_genome".to_owned(),
            read_id_mapping: "explicit_read_index_map".to_owned(),
            source_files: vec!["data.indiv_proba.csv".to_owned()],
        };
        let path = temp_file("assay.json");
        let mut bytes = Vec::new();
        write_assay_metadata_to_writer(&mut bytes, &metadata).unwrap();
        fs::write(&path, bytes).unwrap();
        assert_eq!(read_assay_metadata(&path).unwrap(), metadata);

        fs::write(
            &path,
            r#"{"schema_version":1,"assay_id":"x","caller":"x","caller_version":"x","model_id":"x","chemistry":"x","candidate_rule":"x","source_emission_threshold":null,"source_site_filter":"x","candidate_observations_complete":true,"implicit_skip_policy":"not_applicable","coordinate_source":"x","read_id_mapping":"x","source_files":[],"unexpected":1}"#,
        )
        .unwrap();
        assert!(read_assay_metadata(&path).is_err());
        let _ = fs::remove_file(path);
    }
}
