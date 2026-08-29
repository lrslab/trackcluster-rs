use std::collections::BTreeSet;
use std::ffi::OsString;
use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::Context;
use clap::{Args as ClapArgs, ValueEnum};

use crate::io::modbam::{
    read_dorado_pg_details, read_modbam, DoradoPgProvenance, DoradoPgRecordProvenance,
    InvalidRecordPolicy, MmQuestionMarkPolicy, ModBamImportResult, ModBamOptions, ModBamQc,
};
use crate::modification::{AssayMetadata, ProvenanceStatus, MODIFICATION_SCHEMA_VERSION};

#[derive(Clone, Copy, Debug, ValueEnum)]
enum InvalidRecordPolicyArg {
    Fail,
    Skip,
}

impl From<InvalidRecordPolicyArg> for InvalidRecordPolicy {
    fn from(value: InvalidRecordPolicyArg) -> Self {
        match value {
            InvalidRecordPolicyArg::Fail => Self::Fail,
            InvalidRecordPolicyArg::Skip => Self::Skip,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
enum MmQuestionMarkPolicyArg {
    Unknown,
    BelowEmissionThreshold,
}

impl From<MmQuestionMarkPolicyArg> for MmQuestionMarkPolicy {
    fn from(value: MmQuestionMarkPolicyArg) -> Self {
        match value {
            MmQuestionMarkPolicyArg::Unknown => Self::Unknown,
            MmQuestionMarkPolicyArg::BelowEmissionThreshold => Self::BelowEmissionThreshold,
        }
    }
}

fn parse_probability(value: &str) -> Result<f64, String> {
    let probability = value
        .parse::<f64>()
        .map_err(|_| format!("invalid probability {value:?}"))?;
    if !probability.is_finite() || !(0.0..=1.0).contains(&probability) {
        return Err("probability must be finite and in [0, 1]".to_owned());
    }
    Ok(probability)
}

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Biological sample identifier used to prefix BAM query names.
    #[arg(long)]
    pub sample: String,
    /// Assay compatibility stratum; share only across compatible models and settings.
    #[arg(long)]
    pub assay_id: String,
    /// Primary genome-aligned BAM containing MM, ML, and MN tags.
    #[arg(long)]
    pub bam: PathBuf,
    /// One SAM modification code to import, for example A+a.
    #[arg(long = "mod-code", default_value = "A+a")]
    pub mod_code: String,
    /// Exact Dorado modification model identifier.
    #[arg(long = "model-id")]
    pub model_id: String,
    /// Sequencing chemistry recorded in assay provenance.
    #[arg(long, default_value = "RNA004")]
    pub chemistry: String,
    /// Dorado version recorded in assay provenance.
    #[arg(long = "caller-version", default_value = "unknown")]
    pub caller_version: String,
    /// Candidate universe: all-target-canonical-bases or a centered odd-length IUPAC motif.
    #[arg(long = "candidate-rule", default_value = "all-target-canonical-bases")]
    pub candidate_rule: String,
    /// Dorado --modified-bases-threshold used to create the BAM.
    #[arg(long = "source-emission-threshold", value_parser = parse_probability)]
    pub source_emission_threshold: Option<f64>,
    /// Meaning of candidates omitted from MM groups marked `?`.
    #[arg(
        long = "question-mark-policy",
        value_enum,
        default_value_t = MmQuestionMarkPolicyArg::Unknown
    )]
    question_mark_policy: MmQuestionMarkPolicyArg,
    /// Minimum MAPQ for retained primary alignments.
    #[arg(long = "min-mapq", default_value_t = 0)]
    pub min_mapq: u8,
    /// Fail on an invalid primary record or skip it and mark the universe incomplete.
    #[arg(long = "invalid-record-policy", value_enum, default_value_t = InvalidRecordPolicyArg::Fail)]
    invalid_record_policy: InvalidRecordPolicyArg,
    /// Prefix for observations, assay provenance, and import QC outputs.
    #[arg(short, long)]
    pub out: PathBuf,
}

fn append_suffix(prefix: &Path, suffix: &str) -> PathBuf {
    let mut value: OsString = prefix.as_os_str().to_os_string();
    value.push(suffix);
    PathBuf::from(value)
}

#[derive(Clone, Debug)]
struct ResolvedProvenance {
    caller_version: String,
    source_emission_threshold: Option<f64>,
    status: ProvenanceStatus,
}

fn model_matches(declared: &str, recovered: &str) -> bool {
    fn basename(value: &str) -> &str {
        Path::new(value)
            .file_name()
            .and_then(|name| name.to_str())
            .unwrap_or(value)
    }
    declared == recovered || basename(declared) == basename(recovered)
}

fn resolve_provenance(
    args: &Args,
    provenance: &DoradoPgProvenance,
    records: &[DoradoPgRecordProvenance],
) -> anyhow::Result<ResolvedProvenance> {
    let records_with_models = records
        .iter()
        .filter(|record| !record.model_ids.is_empty())
        .collect::<Vec<_>>();
    let model_matching_records = records_with_models
        .iter()
        .copied()
        .filter(|record| {
            record
                .model_ids
                .iter()
                .any(|recovered| model_matches(&args.model_id, recovered))
        })
        .collect::<Vec<_>>();
    if !records_with_models.is_empty() && model_matching_records.is_empty() {
        anyhow::bail!(
            "declared Dorado model {:?} conflicts with @PG modified-base model argument(s) {:?}",
            args.model_id,
            provenance.model_ids
        );
    }
    let model_verified = !model_matching_records.is_empty();

    let complete_matching_records = model_matching_records
        .iter()
        .copied()
        .filter(|record| {
            record.caller_version.is_some() && record.source_emission_thresholds.len() == 1
        })
        .collect::<Vec<_>>();
    let relevant_records = if !complete_matching_records.is_empty() {
        complete_matching_records
    } else if !records_with_models.is_empty() {
        model_matching_records
    } else {
        records.iter().collect::<Vec<_>>()
    };

    let recovered_versions = relevant_records
        .iter()
        .filter_map(|record| record.caller_version.as_deref())
        .collect::<BTreeSet<_>>();
    if recovered_versions.len() > 1 {
        anyhow::bail!(
            "model-matching Dorado @PG records contain conflicting versions: {:?}",
            recovered_versions
        );
    }
    let mut recovered_thresholds = relevant_records
        .iter()
        .flat_map(|record| record.source_emission_thresholds.iter().copied())
        .collect::<Vec<_>>();
    recovered_thresholds.sort_by(f64::total_cmp);
    recovered_thresholds.dedup_by(|left, right| left.total_cmp(right).is_eq());
    if recovered_thresholds.len() > 1 {
        anyhow::bail!(
            "model-matching Dorado @PG records contain conflicting --modified-bases-threshold values: {:?}",
            recovered_thresholds
        );
    }

    let recovered_version = recovered_versions.first().copied();
    if let Some(recovered) = recovered_version {
        if args.caller_version != "unknown" && args.caller_version != recovered {
            anyhow::bail!(
                "declared Dorado caller version {:?} conflicts with @PG VN {:?}",
                args.caller_version,
                recovered
            );
        }
    }
    let recovered_threshold = recovered_thresholds.first().copied();
    if let (Some(declared), Some(recovered)) = (args.source_emission_threshold, recovered_threshold)
    {
        if declared.total_cmp(&recovered).is_ne() {
            anyhow::bail!(
                "declared Dorado source emission threshold {declared} conflicts with @PG --modified-bases-threshold {recovered}"
            );
        }
    }
    let caller_version = recovered_version
        .map(str::to_owned)
        .unwrap_or_else(|| args.caller_version.clone());
    let source_emission_threshold = recovered_threshold.or(args.source_emission_threshold);
    let fully_verified = relevant_records.iter().any(|record| {
        record.caller_version.is_some() && record.source_emission_thresholds.len() == 1
    }) && recovered_version.is_some()
        && model_verified
        && recovered_threshold.is_some();
    let has_user_declaration =
        args.caller_version != "unknown" || args.source_emission_threshold.is_some();
    let fully_declared =
        has_user_declaration && caller_version != "unknown" && source_emission_threshold.is_some();
    let status = if fully_verified {
        ProvenanceStatus::VerifiedFromPg
    } else if fully_declared {
        ProvenanceStatus::UserDeclared
    } else {
        ProvenanceStatus::Unavailable
    };
    Ok(ResolvedProvenance {
        caller_version,
        source_emission_threshold,
        status,
    })
}

fn optional_values<T: ToString>(values: &[T]) -> String {
    if values.is_empty() {
        "NA".to_owned()
    } else {
        values
            .iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>()
            .join("|")
    }
}

fn write_qc<W: Write>(
    mut writer: W,
    result: &ModBamImportResult,
    provenance: &DoradoPgProvenance,
    resolved: &ResolvedProvenance,
) -> anyhow::Result<()> {
    let qc: &ModBamQc = &result.qc;
    writeln!(writer, "metric\tvalue")?;
    for (metric, value) in [
        ("total_records", qc.total_records),
        ("skipped_unmapped", qc.skipped_unmapped),
        ("skipped_secondary", qc.skipped_secondary),
        ("skipped_supplementary", qc.skipped_supplementary),
        ("skipped_below_mapq", qc.skipped_below_mapq),
        ("duplicate_primary_reads", qc.duplicate_primary_reads),
        (
            "skipped_duplicate_primary_records",
            qc.skipped_duplicate_primary_records,
        ),
        ("retained_records", qc.retained_records),
        ("target_canonical_bases", qc.target_canonical_bases),
        ("target_candidate_bases", qc.target_candidate_bases),
        (
            "records_without_target_group",
            qc.records_without_target_group,
        ),
        (
            "target_groups_skip_flag_omitted",
            qc.target_groups_skip_flag_omitted,
        ),
        (
            "target_groups_low_probability",
            qc.target_groups_low_probability,
        ),
        ("target_groups_unknown", qc.target_groups_unknown),
        (
            "implicit_low_probability_candidates",
            qc.implicit_low_probability_candidates,
        ),
        ("unknown_candidates", qc.unknown_candidates),
        ("ml_values_consumed", qc.ml_values_consumed),
        ("explicit_target_calls", qc.explicit_target_calls),
        ("target_calls_in_insertions", qc.target_calls_in_insertions),
        ("target_calls_in_soft_clips", qc.target_calls_in_soft_clips),
        (
            "implicit_calls_in_insertions",
            qc.implicit_calls_in_insertions,
        ),
        (
            "implicit_calls_in_soft_clips",
            qc.implicit_calls_in_soft_clips,
        ),
        (
            "emitted_explicit_observations",
            qc.emitted_explicit_observations,
        ),
        (
            "emitted_implicit_observations",
            qc.emitted_implicit_observations,
        ),
        (
            "emitted_unknown_observations",
            qc.emitted_unknown_observations,
        ),
        ("emitted_observations", qc.emitted_observations),
    ] {
        writeln!(writer, "{metric}\t{value}")?;
    }
    for (reason, count) in &qc.invalid_records {
        writeln!(writer, "invalid_record_{reason:?}\t{count}")?;
    }
    writeln!(
        writer,
        "candidate_observations_complete\t{}",
        result.semantics.candidate_observations_complete
    )?;
    writeln!(
        writer,
        "implicit_skip_policy\t{}",
        result.semantics.implicit_skip_policy
    )?;
    writeln!(
        writer,
        "source_emission_threshold\t{}",
        result
            .semantics
            .source_emission_threshold
            .map(|value| value.to_string())
            .unwrap_or_else(|| "NA".to_owned())
    )?;
    writeln!(
        writer,
        "mm_question_mark_policy\t{}",
        result.semantics.mm_question_mark_policy.as_str()
    )?;
    writeln!(
        writer,
        "ml_probability_semantics\t{}",
        result.semantics.ml_probability_semantics.as_str()
    )?;
    writeln!(writer, "provenance_status\t{}", resolved.status)?;
    writeln!(
        writer,
        "dorado_pg_records\t{}",
        provenance.dorado_program_records
    )?;
    writeln!(
        writer,
        "dorado_pg_caller_versions\t{}",
        optional_values(&provenance.caller_versions)
    )?;
    writeln!(
        writer,
        "dorado_pg_model_ids\t{}",
        optional_values(&provenance.model_ids)
    )?;
    writeln!(
        writer,
        "dorado_pg_source_emission_thresholds\t{}",
        optional_values(&provenance.source_emission_thresholds)
    )?;
    writeln!(
        writer,
        "dorado_pg_command_lines\t{}",
        optional_values(&provenance.command_lines)
    )?;
    writer.flush().context("flush Dorado/modBAM import QC")?;
    Ok(())
}

fn assay_metadata(
    args: &Args,
    result: &ModBamImportResult,
    resolved: &ResolvedProvenance,
) -> AssayMetadata {
    AssayMetadata {
        schema_version: MODIFICATION_SCHEMA_VERSION,
        assay_id: args.assay_id.clone(),
        caller: "dorado".to_owned(),
        caller_version: resolved.caller_version.clone(),
        model_id: args.model_id.clone(),
        chemistry: args.chemistry.clone(),
        candidate_rule: result.semantics.candidate_rule.clone(),
        source_emission_threshold: result.semantics.source_emission_threshold,
        source_site_filter: format!(
            "primary_genome_aligned_projectable_query_bases;min_mapq={};target={};mm_question_mark_policy={}",
            args.min_mapq,
            result.semantics.target_mod_code,
            result.semantics.mm_question_mark_policy.as_str()
        ),
        candidate_observations_complete: result.semantics.candidate_observations_complete,
        provenance_status: resolved.status,
        implicit_skip_policy: result.semantics.implicit_skip_policy,
        coordinate_source: format!(
            "genome_aligned_bam_cigar_{}",
            result.semantics.ml_probability_semantics.as_str()
        ),
        read_id_mapping: "sample_prefixed_bam_query_name".to_owned(),
        source_files: vec![args.bam.to_string_lossy().into_owned()],
    }
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let observations_path = append_suffix(&args.out, ".observations.tsv");
    let assay_path = append_suffix(&args.out, ".assay.json");
    let qc_path = append_suffix(&args.out, ".import_qc.tsv");
    crate::cli::ensure_distinct_inputs_and_outputs(
        &[("Dorado/modBAM", args.bam.as_path())],
        &[
            (
                "normalized observations output",
                observations_path.as_path(),
            ),
            ("assay metadata output", assay_path.as_path()),
            ("import QC output", qc_path.as_path()),
        ],
    )?;

    let (provenance, provenance_records) = read_dorado_pg_details(&args.bam)?;
    let resolved = resolve_provenance(&args, &provenance, &provenance_records)?;
    let mut options = ModBamOptions::new(&args.assay_id, &args.sample, &args.mod_code);
    options.candidate_rule.clone_from(&args.candidate_rule);
    options.min_mapq = args.min_mapq;
    options.invalid_record_policy = args.invalid_record_policy.into();
    options.source_emission_threshold = resolved.source_emission_threshold;
    options.mm_question_mark_policy = args.question_mark_policy.into();
    let result = read_modbam(&args.bam, &options)?;
    let metadata = assay_metadata(&args, &result, &resolved);
    metadata.validate().map_err(anyhow::Error::msg)?;

    crate::flow::artifact_manifest::atomic_write_with(&observations_path, |writer| {
        crate::io::mod_calls::write_observations_tsv_to_writer(writer, &result.observations)
    })?;
    crate::flow::artifact_manifest::atomic_write_with(&assay_path, |writer| {
        crate::io::mod_calls::write_assay_metadata_to_writer(writer, &metadata)
    })?;
    crate::flow::artifact_manifest::atomic_write_with(&qc_path, |writer| {
        write_qc(writer, &result, &provenance, &resolved)
    })?;

    eprintln!(
        "mod-import-dorado: sample={} assay={} records={} retained={} observations={} complete={}",
        args.sample,
        args.assay_id,
        result.qc.total_records,
        result.qc.retained_records,
        result.qc.emitted_observations,
        result.semantics.candidate_observations_complete
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn args() -> Args {
        Args {
            sample: "S1".to_owned(),
            assay_id: "a1".to_owned(),
            bam: PathBuf::from("calls.bam"),
            mod_code: "A+a".to_owned(),
            model_id: "rna004-test-m6a".to_owned(),
            chemistry: "RNA004".to_owned(),
            caller_version: "unknown".to_owned(),
            candidate_rule: "all-target-canonical-bases".to_owned(),
            source_emission_threshold: None,
            question_mark_policy: MmQuestionMarkPolicyArg::Unknown,
            min_mapq: 0,
            invalid_record_policy: InvalidRecordPolicyArg::Fail,
            out: PathBuf::from("out"),
        }
    }

    fn record(
        caller_version: Option<&str>,
        model_ids: &[&str],
        thresholds: &[f64],
    ) -> DoradoPgRecordProvenance {
        DoradoPgRecordProvenance {
            caller_version: caller_version.map(str::to_owned),
            model_ids: model_ids.iter().map(|value| (*value).to_owned()).collect(),
            source_emission_thresholds: thresholds.to_vec(),
            command_line: None,
        }
    }

    #[test]
    fn probability_parser_is_strict() {
        assert_eq!(parse_probability("0").unwrap(), 0.0);
        assert_eq!(parse_probability("1").unwrap(), 1.0);
        assert!(parse_probability("NaN").is_err());
        assert!(parse_probability("1.1").is_err());
    }

    #[test]
    fn provenance_is_not_verified_by_combining_different_pg_records() {
        let records = vec![
            record(Some("0.9.1"), &[], &[]),
            record(None, &["rna004-test-m6a"], &[0.05]),
        ];
        let provenance = DoradoPgProvenance {
            dorado_program_records: 2,
            caller_versions: vec!["0.9.1".to_owned()],
            model_ids: vec!["rna004-test-m6a".to_owned()],
            source_emission_thresholds: vec![0.05],
            command_lines: Vec::new(),
        };

        let resolved = resolve_provenance(&args(), &provenance, &records).unwrap();
        assert_eq!(resolved.caller_version, "unknown");
        assert_eq!(resolved.source_emission_threshold, Some(0.05));
        assert_eq!(resolved.status, ProvenanceStatus::Unavailable);
    }

    #[test]
    fn complete_model_record_is_not_conflicted_by_unrelated_dorado_stage() {
        let records = vec![
            record(Some("0.9.1"), &["rna004-test-m6a"], &[0.05]),
            record(Some("1.2.0"), &[], &[]),
        ];
        let provenance = DoradoPgProvenance {
            dorado_program_records: 2,
            caller_versions: vec!["0.9.1".to_owned(), "1.2.0".to_owned()],
            model_ids: vec!["rna004-test-m6a".to_owned()],
            source_emission_thresholds: vec![0.05],
            command_lines: Vec::new(),
        };

        let resolved = resolve_provenance(&args(), &provenance, &records).unwrap();
        assert_eq!(resolved.caller_version, "0.9.1");
        assert_eq!(resolved.source_emission_threshold, Some(0.05));
        assert_eq!(resolved.status, ProvenanceStatus::VerifiedFromPg);
    }
}
