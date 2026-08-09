use std::ffi::OsString;
use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::Context;
use clap::Args as ClapArgs;

use crate::io::gff::{AnnotationFormat, GffToBiggOptions};
use crate::io::m6anet::{import_m6anet_with_site_probability, M6anetImportOptions, M6anetImportQc};

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Biological sample identifier used to prefix normalized read IDs.
    #[arg(long)]
    pub sample: String,
    /// Assay compatibility stratum; use the same value only for compatible runs.
    #[arg(long)]
    pub assay_id: String,
    /// m6Anet data.indiv_proba.csv or data.indiv_proba.csv.gz.
    #[arg(long = "indiv")]
    pub indiv: PathBuf,
    /// Optional m6Anet data.info used to audit the retained candidate universe.
    #[arg(long = "data-info")]
    pub data_info: Option<PathBuf>,
    /// Optional data.site_proba.csv used only for retained-site and mod-ratio QC.
    #[arg(long = "site-proba")]
    pub site_proba: Option<PathBuf>,
    /// Explicit TSV mapping with columns read_index and read_id.
    #[arg(long = "read-map")]
    pub read_map: PathBuf,
    /// Reference GTF/GFF whose transcript IDs exactly match m6Anet output.
    #[arg(long = "reference")]
    pub reference: PathBuf,
    /// Exact m6Anet model identifier.
    #[arg(long = "model-id")]
    pub model_id: String,
    /// m6Anet version recorded in assay provenance.
    #[arg(long = "caller-version", default_value = "unknown")]
    pub caller_version: String,
    /// Centered odd-length IUPAC candidate motif with A at the target position.
    #[arg(long = "candidate-rule", default_value = "DRACH")]
    pub candidate_rule: String,
    /// Expected data.info minimum read filter; enables exact retained-site audit.
    #[arg(long = "min-reads")]
    pub min_reads: Option<u64>,
    /// Read threshold used to reproduce source mod_ratio; known model presets are automatic.
    #[arg(long = "read-probability-threshold")]
    pub read_probability_threshold: Option<f64>,
    /// Reference annotation attribute syntax.
    #[arg(long = "input-format", value_enum, default_value_t = AnnotationFormat::Auto)]
    pub input_format: AnnotationFormat,
    /// GFF3 gene-feature attribute written as the projection catalog gene ID.
    #[arg(long = "gene-key", default_value = "ID")]
    pub gene_key: String,
    /// Prefix for observations, assay provenance, and import QC outputs.
    #[arg(short, long)]
    pub out: PathBuf,
}

fn append_suffix(prefix: &Path, suffix: &str) -> PathBuf {
    let mut value: OsString = prefix.as_os_str().to_os_string();
    value.push(suffix);
    PathBuf::from(value)
}

fn write_qc<W: Write>(mut writer: W, qc: &M6anetImportQc) -> anyhow::Result<()> {
    writeln!(writer, "metric\tvalue")?;
    writeln!(writer, "input_rows\t{}", qc.input_rows)?;
    writeln!(writer, "unique_observations\t{}", qc.unique_observations)?;
    writeln!(writer, "duplicate_exact\t{}", qc.duplicate_exact)?;
    writeln!(writer, "read_map_entries\t{}", qc.read_map_entries)?;
    writeln!(
        writer,
        "read_map_entries_used\t{}",
        qc.read_map_entries_used
    )?;
    writeln!(writer, "source_transcripts\t{}", qc.source_transcripts)?;
    writeln!(
        writer,
        "projection_transcripts_loaded\t{}",
        qc.projection_transcripts_loaded
    )?;
    writeln!(writer, "source_sites\t{}", qc.source_sites)?;
    if let Some(info) = &qc.data_info {
        writeln!(writer, "data_info_sites\t{}", info.sites)?;
        writeln!(writer, "data_info_retained_sites\t{}", info.retained_sites)?;
        writeln!(writer, "data_info_filtered_sites\t{}", info.filtered_sites)?;
        writeln!(writer, "data_info_total_reads\t{}", info.total_reads)?;
        writeln!(writer, "data_info_retained_reads\t{}", info.retained_reads)?;
        writeln!(
            writer,
            "data_info_minimum_reads\t{}",
            info.minimum_reads
                .map(|value| value.to_string())
                .unwrap_or_else(|| "NA".to_owned())
        )?;
    } else {
        writeln!(writer, "data_info_sites\tNA")?;
        writeln!(writer, "data_info_retained_sites\tNA")?;
        writeln!(writer, "data_info_filtered_sites\tNA")?;
        writeln!(writer, "data_info_total_reads\tNA")?;
        writeln!(writer, "data_info_retained_reads\tNA")?;
        writeln!(writer, "data_info_minimum_reads\tNA")?;
    }
    if let Some(site) = &qc.site_probability {
        writeln!(writer, "site_probability_sites\t{}", site.sites)?;
        writeln!(writer, "site_probability_total_reads\t{}", site.total_reads)?;
        writeln!(
            writer,
            "site_probability_sites_at_or_above_threshold\t{}",
            site.sites_at_or_above_probability_threshold
        )?;
        writeln!(
            writer,
            "site_probability_threshold\t{}",
            site.site_probability_threshold
        )?;
        writeln!(
            writer,
            "read_probability_threshold\t{}",
            site.read_probability_threshold
                .map(|value| value.to_string())
                .unwrap_or_else(|| "NA".to_owned())
        )?;
    } else {
        writeln!(writer, "site_probability_sites\tNA")?;
        writeln!(writer, "site_probability_total_reads\tNA")?;
        writeln!(writer, "site_probability_sites_at_or_above_threshold\tNA")?;
        writeln!(writer, "site_probability_threshold\tNA")?;
        writeln!(writer, "read_probability_threshold\tNA")?;
    }
    writer.flush().context("flush m6Anet import QC")?;
    Ok(())
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let observations_path = append_suffix(&args.out, ".observations.tsv");
    let assay_path = append_suffix(&args.out, ".assay.json");
    let qc_path = append_suffix(&args.out, ".import_qc.tsv");
    let mut inputs = vec![
        ("m6Anet individual probabilities", args.indiv.as_path()),
        ("m6Anet read map", args.read_map.as_path()),
        ("projection reference", args.reference.as_path()),
    ];
    if let Some(path) = args.data_info.as_deref() {
        inputs.push(("m6Anet data.info", path));
    }
    if let Some(path) = args.site_proba.as_deref() {
        inputs.push(("m6Anet data.site_proba", path));
    }
    crate::cli::ensure_distinct_inputs_and_outputs(
        &inputs,
        &[
            (
                "normalized observations output",
                observations_path.as_path(),
            ),
            ("assay metadata output", assay_path.as_path()),
            ("import QC output", qc_path.as_path()),
        ],
    )?;

    let mut options = M6anetImportOptions::new(&args.sample, &args.assay_id, &args.model_id);
    options.caller_version = args.caller_version;
    options.candidate_rule = args.candidate_rule;
    options.minimum_reads = args.min_reads;
    if let Some(threshold) = args.read_probability_threshold {
        options.read_probability_threshold = Some(threshold);
    }
    options.annotation_options = GffToBiggOptions {
        format: args.input_format,
        gene_key: args.gene_key,
    };
    let result = import_m6anet_with_site_probability(
        &args.indiv,
        &args.read_map,
        &args.reference,
        args.data_info.as_deref(),
        args.site_proba.as_deref(),
        &options,
    )?;

    crate::flow::artifact_manifest::atomic_write_with(&observations_path, |writer| {
        crate::io::mod_calls::write_observations_tsv_to_writer(writer, &result.observations)
    })?;
    crate::flow::artifact_manifest::atomic_write_with(&assay_path, |writer| {
        crate::io::mod_calls::write_assay_metadata_to_writer(writer, &result.metadata)
    })?;
    crate::flow::artifact_manifest::atomic_write_with(&qc_path, |writer| {
        write_qc(writer, &result.qc)
    })?;

    eprintln!(
        "mod-import-m6anet: sample={} assay={} input_rows={} observations={} sites={}",
        args.sample,
        args.assay_id,
        result.qc.input_rows,
        result.qc.unique_observations,
        result.qc.source_sites
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn qc_writer_reports_absent_data_info_as_na() {
        let mut output = Vec::new();
        write_qc(&mut output, &M6anetImportQc::default()).unwrap();
        let output = String::from_utf8(output).unwrap();
        assert!(output.starts_with("metric\tvalue\ninput_rows\t0\n"));
        assert!(output.contains("data_info_sites\tNA\n"));
    }
}
