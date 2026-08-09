use std::collections::{BTreeMap, BTreeSet};
use std::ffi::OsString;
use std::path::{Path, PathBuf};

use anyhow::Context;
use clap::Args as ClapArgs;

use crate::modification::aggregate::{
    aggregate_modifications, write_isoform_mod_design_tsv, write_isoform_mod_sites_tsv,
    write_join_qc_tsv, write_site_join_qc_tsv, AggregateOptions, ModSampleInput,
};

#[derive(Clone, Debug)]
pub(crate) struct AnalysisThreshold {
    pub(crate) assay_id: String,
    pub(crate) value: f64,
}

pub(crate) fn parse_analysis_threshold(value: &str) -> Result<AnalysisThreshold, String> {
    let (assay_id, threshold) = value
        .split_once('=')
        .ok_or_else(|| "expected ASSAY_ID=PROBABILITY".to_owned())?;
    if assay_id.trim().is_empty() || assay_id.chars().any(char::is_control) {
        return Err("assay id must not be empty or contain control characters".to_owned());
    }
    let threshold = threshold
        .parse::<f64>()
        .map_err(|_| format!("invalid probability {threshold:?}"))?;
    if !threshold.is_finite() || !(0.0..=1.0).contains(&threshold) {
        return Err("probability must be finite and in [0, 1]".to_owned());
    }
    Ok(AnalysisThreshold {
        assay_id: assay_id.to_owned(),
        value: threshold,
    })
}

pub(crate) fn parse_rate(value: &str) -> Result<f64, String> {
    let rate = value
        .parse::<f64>()
        .map_err(|_| format!("invalid rate {value:?}"))?;
    if !rate.is_finite() || !(0.0..=1.0).contains(&rate) {
        return Err("rate must be finite and in [0, 1]".to_owned());
    }
    Ok(rate)
}

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// TrackCluster samples manifest containing sample, group, and reads columns.
    #[arg(long)]
    pub manifest: PathBuf,
    /// Final pooled isoform BED/bigGenePred catalog.
    #[arg(long)]
    pub isoforms: PathBuf,
    /// Final unique read-to-isoform TSV.
    #[arg(long = "read-to-isoform")]
    pub read_to_isoform: PathBuf,
    /// Modification samples manifest.
    #[arg(long = "mod-manifest")]
    pub mod_manifest: PathBuf,
    /// Indexed genomic FASTA used to reject canonical-base mismatches.
    #[arg(long = "reference-fasta")]
    pub reference_fasta: Option<PathBuf>,
    /// Per-assay hard-call threshold as ASSAY_ID=PROBABILITY; repeat for each assay.
    #[arg(long = "analysis-threshold", required = true, value_parser = parse_analysis_threshold)]
    analysis_thresholds: Vec<AnalysisThreshold>,
    /// Minimum callable molecules for an eligible sample/isoform/site row.
    #[arg(long, default_value_t = 1, value_parser = clap::value_parser!(u64).range(1..))]
    pub min_callable: u64,
    /// Minimum distinct-read join rate for each sample/assay.
    #[arg(long, default_value = "0.9", value_parser = parse_rate)]
    pub min_read_join_rate: f64,
    /// Emit low-join rows as ineligible instead of failing.
    #[arg(long)]
    pub allow_low_join: bool,
    /// Output prefix for the four modification aggregate TSV files.
    #[arg(short, long)]
    pub out: PathBuf,
}

fn append_suffix(prefix: &Path, suffix: &str) -> PathBuf {
    let mut value: OsString = prefix.as_os_str().to_os_string();
    value.push(suffix);
    PathBuf::from(value)
}

fn ensure_nested_inputs_are_not_outputs(
    rows: &[crate::io::mod_manifest::ModSampleRow],
    outputs: &[(&str, &Path)],
) -> anyhow::Result<()> {
    let mut labeled_inputs = Vec::new();
    for row in rows {
        labeled_inputs.push((
            format!(
                "modification observations for sample {:?}, assay {:?}",
                row.sample, row.assay_id
            ),
            row.observations.as_path(),
        ));
        labeled_inputs.push((
            format!(
                "assay metadata for sample {:?}, assay {:?}",
                row.sample, row.assay_id
            ),
            row.assay_metadata.as_path(),
        ));
        if let Some(path) = row.coverage_bam.as_deref() {
            labeled_inputs.push((
                format!(
                    "coverage BAM for sample {:?}, assay {:?}",
                    row.sample, row.assay_id
                ),
                path,
            ));
        }
    }
    let inputs = labeled_inputs
        .iter()
        .map(|(label, path)| (label.as_str(), *path))
        .collect::<Vec<_>>();
    crate::cli::ensure_distinct_inputs_and_outputs(&inputs, outputs)
}

pub(crate) fn collect_thresholds(
    analysis_thresholds: Vec<AnalysisThreshold>,
) -> anyhow::Result<BTreeMap<String, f64>> {
    let mut thresholds = BTreeMap::new();
    for threshold in analysis_thresholds {
        if thresholds
            .insert(threshold.assay_id.clone(), threshold.value)
            .is_some()
        {
            anyhow::bail!(
                "analysis threshold supplied more than once for assay {:?}",
                threshold.assay_id
            );
        }
    }
    Ok(thresholds)
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn run_with_paths(
    manifest: &Path,
    isoforms: &Path,
    read_to_isoform: &Path,
    mod_manifest: &Path,
    reference_fasta: Option<&Path>,
    thresholds: BTreeMap<String, f64>,
    min_callable: u64,
    min_read_join_rate: f64,
    allow_low_join: bool,
    out: &Path,
) -> anyhow::Result<crate::modification::aggregate::AggregateResult> {
    let join_qc_path = append_suffix(out, ".mod_join_qc.tsv");
    let site_join_qc_path = append_suffix(out, ".mod_site_join_qc.tsv");
    let sites_path = append_suffix(out, ".isoform_mod_sites.tsv");
    let design_path = append_suffix(out, ".isoform_mod_design.tsv");
    let mut inputs = vec![
        ("sample manifest", manifest),
        ("isoform catalog", isoforms),
        ("read-to-isoform", read_to_isoform),
        ("mod manifest", mod_manifest),
    ];
    if let Some(path) = reference_fasta {
        inputs.push(("reference FASTA", path));
    }
    crate::cli::ensure_distinct_inputs_and_outputs(
        &inputs,
        &[
            ("join QC output", join_qc_path.as_path()),
            ("site join QC output", site_join_qc_path.as_path()),
            ("isoform modification sites output", sites_path.as_path()),
            ("isoform modification design output", design_path.as_path()),
        ],
    )?;

    let samples = crate::io::manifest::read_manifest_tsv(manifest)?;
    let isoform_records = crate::io::bed::read_bed12(isoforms)
        .with_context(|| format!("open isoforms {isoforms:?}"))?
        .collect::<Result<Vec<_>, _>>()
        .with_context(|| format!("parse isoforms {isoforms:?}"))?;
    let read_to_isoform_records = crate::count::read_read_to_isoform_tsv(read_to_isoform)
        .with_context(|| format!("read unique mapping {read_to_isoform:?}"))?;
    let mod_rows = crate::io::mod_manifest::read_mod_manifest_tsv(mod_manifest)?;
    ensure_nested_inputs_are_not_outputs(
        &mod_rows,
        &[
            ("join QC output", join_qc_path.as_path()),
            ("site join QC output", site_join_qc_path.as_path()),
            ("isoform modification sites output", sites_path.as_path()),
            ("isoform modification design output", design_path.as_path()),
        ],
    )?;

    let mut mod_inputs = Vec::with_capacity(mod_rows.len());
    for row in mod_rows {
        let metadata = crate::io::mod_calls::read_assay_metadata(&row.assay_metadata)?;
        let observations = crate::io::mod_calls::read_observations_tsv(&row.observations)?;
        let coverage = row
            .coverage_bam
            .as_deref()
            .map(|path| crate::io::coverage::read_primary_bam_coverage(path, &row.sample))
            .transpose()?;
        mod_inputs.push(ModSampleInput {
            sample: row.sample,
            assay_id: row.assay_id,
            metadata,
            observations,
            coverage,
        });
    }

    let reference_bases = if let Some(path) = reference_fasta {
        let sites = mod_inputs
            .iter()
            .flat_map(|input| {
                input
                    .observations
                    .observations
                    .iter()
                    .map(|observation| observation.key.site.clone())
            })
            .collect::<BTreeSet<_>>();
        let mut fasta = crate::io::fasta::IndexedFasta::open(path)?;
        let mut bases = BTreeMap::new();
        for site in sites {
            let base = fasta.oriented_base(&site.chrom, site.pos0, site.strand)?;
            bases.insert(site, base);
        }
        Some(bases)
    } else {
        None
    };

    let result = aggregate_modifications(
        &samples,
        &isoform_records,
        &read_to_isoform_records,
        &mod_inputs,
        &AggregateOptions {
            analysis_thresholds: thresholds,
            min_callable,
            min_read_join_rate,
            allow_low_join,
            reference_bases,
        },
    )?;

    crate::flow::artifact_manifest::atomic_write_with(&join_qc_path, |writer| {
        write_join_qc_tsv(writer, &result.join_qc)
    })?;
    crate::flow::artifact_manifest::atomic_write_with(&site_join_qc_path, |writer| {
        write_site_join_qc_tsv(writer, &result.site_join_qc)
    })?;
    crate::flow::artifact_manifest::atomic_write_with(&sites_path, |writer| {
        write_isoform_mod_sites_tsv(writer, &result.sites)
    })?;
    crate::flow::artifact_manifest::atomic_write_with(&design_path, |writer| {
        write_isoform_mod_design_tsv(writer, &result.design)
    })?;

    Ok(result)
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let thresholds = collect_thresholds(args.analysis_thresholds)?;
    let result = run_with_paths(
        &args.manifest,
        &args.isoforms,
        &args.read_to_isoform,
        &args.mod_manifest,
        args.reference_fasta.as_deref(),
        thresholds,
        args.min_callable,
        args.min_read_join_rate,
        args.allow_low_join,
        &args.out,
    )?;

    eprintln!(
        "mod-aggregate: assays={} join_qc_rows={} site_join_qc_rows={} site_rows={} design_rows={}",
        result
            .join_qc
            .iter()
            .map(|row| row.assay_id.as_str())
            .collect::<std::collections::BTreeSet<_>>()
            .len(),
        result.join_qc.len(),
        result.site_join_qc.len(),
        result.sites.len(),
        result.design.len()
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::time::{SystemTime, UNIX_EPOCH};

    use super::*;

    #[test]
    fn threshold_parser_rejects_missing_duplicate_separator_and_nonfinite_values() {
        assert!(parse_analysis_threshold("a=0").is_ok());
        assert!(parse_analysis_threshold("a=1").is_ok());
        assert!(parse_analysis_threshold("a").is_err());
        assert!(parse_analysis_threshold("=0.5").is_err());
        assert!(parse_analysis_threshold("a=NaN").is_err());
        assert!(parse_analysis_threshold("a=1.1").is_err());
    }

    #[test]
    fn nested_manifest_inputs_cannot_alias_generated_outputs() {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let root = std::env::temp_dir().join(format!(
            "trackcluster-mod-aggregate-alias-{}-{nonce}",
            std::process::id()
        ));
        fs::create_dir(&root).unwrap();
        let observations = root.join("result.isoform_mod_design.tsv");
        let metadata = root.join("assay.json");
        fs::write(&observations, "input\n").unwrap();
        fs::write(&metadata, "{}\n").unwrap();
        let rows = vec![crate::io::mod_manifest::ModSampleRow {
            sample: "S1".to_owned(),
            assay_id: "a1".to_owned(),
            observations: observations.clone(),
            assay_metadata: metadata,
            coverage_bam: None,
        }];
        let error = ensure_nested_inputs_are_not_outputs(
            &rows,
            &[("design output", observations.as_path())],
        )
        .unwrap_err()
        .to_string();
        assert!(
            error.contains("same path") || error.contains("same file") || error.contains("aliases"),
            "{error}"
        );
        assert_eq!(fs::read_to_string(&observations).unwrap(), "input\n");
        fs::remove_dir_all(root).unwrap();
    }
}
