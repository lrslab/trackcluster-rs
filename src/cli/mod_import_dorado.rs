use std::ffi::OsString;
use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::Context;
use clap::{Args as ClapArgs, ValueEnum};

use crate::io::modbam::{
    read_modbam, InvalidRecordPolicy, MmQuestionMarkPolicy, ModBamImportResult, ModBamOptions,
    ModBamQc,
};
use crate::modification::{AssayMetadata, MODIFICATION_SCHEMA_VERSION};

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

fn write_qc<W: Write>(mut writer: W, result: &ModBamImportResult) -> anyhow::Result<()> {
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
    writer.flush().context("flush Dorado/modBAM import QC")?;
    Ok(())
}

fn assay_metadata(args: &Args, result: &ModBamImportResult) -> AssayMetadata {
    AssayMetadata {
        schema_version: MODIFICATION_SCHEMA_VERSION,
        assay_id: args.assay_id.clone(),
        caller: "dorado".to_owned(),
        caller_version: args.caller_version.clone(),
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

    let mut options = ModBamOptions::new(&args.assay_id, &args.sample, &args.mod_code);
    options.candidate_rule.clone_from(&args.candidate_rule);
    options.min_mapq = args.min_mapq;
    options.invalid_record_policy = args.invalid_record_policy.into();
    options.source_emission_threshold = args.source_emission_threshold;
    options.mm_question_mark_policy = args.question_mark_policy.into();
    let result = read_modbam(&args.bam, &options)?;
    let metadata = assay_metadata(&args, &result);
    metadata.validate().map_err(anyhow::Error::msg)?;

    crate::flow::artifact_manifest::atomic_write_with(&observations_path, |writer| {
        crate::io::mod_calls::write_observations_tsv_to_writer(writer, &result.observations)
    })?;
    crate::flow::artifact_manifest::atomic_write_with(&assay_path, |writer| {
        crate::io::mod_calls::write_assay_metadata_to_writer(writer, &metadata)
    })?;
    crate::flow::artifact_manifest::atomic_write_with(&qc_path, |writer| {
        write_qc(writer, &result)
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

    #[test]
    fn probability_parser_is_strict() {
        assert_eq!(parse_probability("0").unwrap(), 0.0);
        assert_eq!(parse_probability("1").unwrap(), 1.0);
        assert!(parse_probability("NaN").is_err());
        assert!(parse_probability("1.1").is_err());
    }
}
