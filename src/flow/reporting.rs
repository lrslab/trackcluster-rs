//! Atomic publication of batch summaries, error logs, and downsampling reports.

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::Context;

use crate::flow::artifact_manifest::{
    atomic_write_with, read_run_manifest, ToolIdentity, MANIFEST_FILE_NAME,
};
use crate::flow::executor::ExecutionReport;
use crate::flow::full::{write_downsample_records_to_writer, BatchRunOptions};
use crate::flow::path_key::{ensure_destination_within, gene_dir_path, GeneId};

/// Immutable context needed to publish one batch report set.
pub(super) struct BatchReportContext<'a> {
    pub(super) args: &'a BatchRunOptions,
    pub(super) batch_file_prefix: &'a str,
    pub(super) gene_path_map: &'a Path,
    pub(super) total: usize,
}

/// Paths and first error produced while publishing a batch report set.
pub(super) struct PublishedBatchReport {
    pub(super) summary_path: PathBuf,
    pub(super) error_path: PathBuf,
    pub(super) downsample_path: PathBuf,
    pub(super) first_error: Option<String>,
}

fn effective_options_json(args: &BatchRunOptions, effective_threads: usize) -> serde_json::Value {
    serde_json::json!({
        "schema_version": 2,
        "cluster_mode": args.cluster_mode.as_str(),
        "prepare": {
            "fraction_read": args.prepare.fraction_read,
            "fraction_ref": args.prepare.fraction_ref,
        },
        "clustering": {
            "sw_score": args.clustering.sw_score,
            "batch_size": args.clustering.batch_size,
            "batch_rounds": args.clustering.batch_rounds,
            "name2_mode": args.clustering.name2_mode.as_str(),
            "junction": {
                "platform_preset": args.clustering.junction.platform_preset.as_str(),
                "correction_offset": args.clustering.junction.correction.offset,
                "correction_min_support": args.clustering.junction.correction.min_support,
                "sl_partial_five_prime_offset": args.clustering.junction.sl.partial_five_prime_end_offset,
                "sl_same_junction_five_prime_offset": args.clustering.junction.sl.same_junction_five_prime_end_offset,
                "sl_five_prime_cluster_offset": args.clustering.junction.sl.five_prime_cluster_offset,
                "sl_five_prime_min_support": args.clustering.junction.sl.min_five_prime_cluster_support,
                "same_junction_three_prime_offset": args.clustering.junction.three_prime.same_junction_three_prime_end_offset,
                "three_prime_cluster_offset": args.clustering.junction.three_prime.three_prime_cluster_offset,
                "three_prime_min_support": args.clustering.junction.three_prime.min_three_prime_cluster_support,
            },
            "overlap": {
                "cutoff1": args.clustering.overlap.cutoff1,
                "cutoff2": args.clustering.overlap.cutoff2,
                "intron_weight": args.clustering.overlap.intron_weight,
            },
        },
        "counting": {
            "assignment_mode": args.counting.assignment_mode.to_string(),
            "unique_assignment_junction_offset": args.counting.unique_assignment.junction_offset,
        },
        "runtime": {
            "requested_threads": args.runtime.threads,
            "effective_threads": effective_threads,
            "force": args.runtime.force,
            "gene_error_policy": args.runtime.gene_error_policy.to_string(),
            "invalid_read_policy": args.runtime.invalid_read_policy.to_string(),
            "progress_every": args.runtime.progress_every,
            "heartbeat_seconds": args.runtime.heartbeat_seconds,
            "heartbeat_top": args.runtime.heartbeat_top,
        },
        "downsample": {
            "genes": args.downsample.genes,
            "max_reads_per_gene": args.downsample.max_reads_per_gene,
            "seed": args.downsample.seed,
        },
    })
}

fn write_summary<W: Write>(
    summary: &mut W,
    context: &BatchReportContext<'_>,
    execution: &ExecutionReport,
) -> anyhow::Result<()> {
    let args = context.args;
    let invocation: Vec<String> = std::env::args_os()
        .map(|argument| argument.to_string_lossy().into_owned())
        .collect();
    let tool = ToolIdentity::current();
    writeln!(
        summary,
        "invocation_json\t{}",
        serde_json::to_string(&invocation)?
    )?;
    writeln!(
        summary,
        "effective_options_json\t{}",
        serde_json::to_string(&effective_options_json(args, execution.worker_count))?
    )?;
    writeln!(summary, "executable_version\t{}", tool.package_version)?;
    writeln!(summary, "git_commit\t{}", tool.git_commit)?;
    writeln!(summary, "source_fingerprint\t{}", tool.source_fingerprint)?;
    writeln!(summary, "input_root\t{:?}", args.input_root)?;
    writeln!(summary, "gene_list\t{:?}", args.gene_list)?;
    writeln!(summary, "output_root\t{:?}", args.output_root)?;
    writeln!(summary, "gene_path_map\t{:?}", context.gene_path_map)?;
    writeln!(summary, "cluster_mode\t{}", args.cluster_mode)?;
    writeln!(summary, "threads\t{}", args.runtime.threads)?;
    writeln!(summary, "effective_threads\t{}", execution.worker_count)?;
    writeln!(
        summary,
        "gene_error_policy\t{}",
        args.runtime.gene_error_policy
    )?;
    writeln!(
        summary,
        "invalid_read_policy\t{}",
        args.runtime.invalid_read_policy
    )?;
    writeln!(
        summary,
        "prepare_fraction_read\t{}",
        args.prepare.fraction_read
    )?;
    writeln!(
        summary,
        "prepare_fraction_ref\t{}",
        args.prepare.fraction_ref
    )?;
    writeln!(summary, "sw_score\t{}", args.clustering.sw_score)?;
    writeln!(summary, "batch_size\t{}", args.clustering.batch_size)?;
    writeln!(summary, "batch_rounds\t{}", args.clustering.batch_rounds)?;
    writeln!(summary, "name2_mode\t{}", args.clustering.name2_mode)?;
    writeln!(
        summary,
        "platform_preset\t{}",
        args.clustering.junction.platform_preset
    )?;
    writeln!(
        summary,
        "junction_correction_offset\t{}",
        args.clustering.junction.correction.offset
    )?;
    writeln!(
        summary,
        "junction_correction_min_support\t{}",
        args.clustering.junction.correction.min_support
    )?;
    writeln!(
        summary,
        "sl_partial_5prime_offset\t{}",
        args.clustering.junction.sl.partial_five_prime_end_offset
    )?;
    writeln!(
        summary,
        "sl_same_junction_5prime_offset\t{}",
        args.clustering
            .junction
            .sl
            .same_junction_five_prime_end_offset
    )?;
    writeln!(
        summary,
        "sl_5prime_cluster_offset\t{}",
        args.clustering.junction.sl.five_prime_cluster_offset
    )?;
    writeln!(
        summary,
        "sl_5prime_min_support\t{}",
        args.clustering.junction.sl.min_five_prime_cluster_support
    )?;
    writeln!(
        summary,
        "same_junction_3prime_offset\t{}",
        args.clustering
            .junction
            .three_prime
            .same_junction_three_prime_end_offset
    )?;
    writeln!(
        summary,
        "3prime_cluster_offset\t{}",
        args.clustering
            .junction
            .three_prime
            .three_prime_cluster_offset
    )?;
    writeln!(
        summary,
        "3prime_min_support\t{}",
        args.clustering
            .junction
            .three_prime
            .min_three_prime_cluster_support
    )?;
    writeln!(
        summary,
        "overlap_cutoff1\t{}",
        args.clustering.overlap.cutoff1
    )?;
    writeln!(
        summary,
        "overlap_cutoff2\t{}",
        args.clustering.overlap.cutoff2
    )?;
    writeln!(
        summary,
        "overlap_intron_weight\t{}",
        args.clustering.overlap.intron_weight
    )?;
    writeln!(
        summary,
        "assignment_mode\t{}",
        args.counting.assignment_mode
    )?;
    writeln!(
        summary,
        "unique_assignment_junction_offset\t{}",
        args.counting.unique_assignment.junction_offset
    )?;
    writeln!(summary, "force\t{}", args.runtime.force)?;
    writeln!(summary, "progress_every\t{}", args.runtime.progress_every)?;
    writeln!(
        summary,
        "heartbeat_seconds\t{}",
        args.runtime.heartbeat_seconds
    )?;
    writeln!(summary, "heartbeat_top\t{}", args.runtime.heartbeat_top)?;
    writeln!(
        summary,
        "max_reads_per_gene\t{}",
        args.downsample.max_reads_per_gene
    )?;
    writeln!(summary, "downsample_seed\t{}", args.downsample.seed)?;
    if args.downsample.genes.is_empty() {
        writeln!(summary, "downsample_genes\t[]")?;
    } else {
        writeln!(
            summary,
            "downsample_genes\t{}",
            args.downsample.genes.join(",")
        )?;
    }
    writeln!(summary, "total_genes\t{}", context.total)?;
    let status = if execution.errors == 0 {
        "complete"
    } else if args.runtime.gene_error_policy == crate::flow::config::GeneErrorPolicy::Continue
        && !execution.mergeable_genes.is_empty()
        && execution.infrastructure_error_count() == 0
    {
        "partial"
    } else {
        "failed"
    };
    writeln!(summary, "status\t{status}")?;
    writeln!(summary, "processed\t{}", execution.processed)?;
    writeln!(summary, "skipped\t{}", execution.skipped)?;
    writeln!(
        summary,
        "skipped_completed_outputs\t{}",
        execution.skipped_completed_outputs
    )?;
    writeln!(
        summary,
        "skipped_empty_reads\t{}",
        execution.skipped_empty_reads
    )?;
    writeln!(
        summary,
        "skipped_no_usable_reads\t{}",
        execution.skipped_no_usable_reads
    )?;
    writeln!(
        summary,
        "all_reads_rejected_genes\t{}",
        execution.skipped_no_usable_reads
    )?;
    writeln!(
        summary,
        "prepare_rejected_read_tracks\t{}",
        execution.prepare_rejected_read_tracks
    )?;
    writeln!(
        summary,
        "per_gene_rejected_read_tracks\t{}",
        execution.rejected_read_tracks
    )?;
    writeln!(
        summary,
        "rejected_read_tracks\t{}",
        execution
            .prepare_rejected_read_tracks
            .saturating_add(execution.rejected_read_tracks)
    )?;
    writeln!(
        summary,
        "genes_with_rejected_reads\t{}",
        execution.genes_with_rejected_reads
    )?;
    writeln!(summary, "errors\t{}", execution.errors)?;
    writeln!(
        summary,
        "failed_missing_inputs\t{}",
        execution.failed_missing_inputs
    )?;
    writeln!(
        summary,
        "failed_processing\t{}",
        execution.failed_processing
    )?;
    writeln!(summary, "failed_panics\t{}", execution.failed_panics)?;
    writeln!(
        summary,
        "mergeable_genes\t{}",
        execution.mergeable_genes.len()
    )?;
    writeln!(
        summary,
        "excluded_failed_genes\t{}",
        execution.failed_gene_count()
    )?;
    writeln!(
        summary,
        "infrastructure_errors\t{}",
        execution.infrastructure_error_count()
    )?;
    writeln!(
        summary,
        "resume_reused\t{}",
        execution
            .resume_decisions
            .iter()
            .filter(|decision| decision.action == "reuse")
            .count()
    )?;
    writeln!(
        summary,
        "resume_rebuilt\t{}",
        execution
            .resume_decisions
            .iter()
            .filter(|decision| decision.action == "rebuild")
            .count()
    )?;
    for decision in &execution.resume_decisions {
        writeln!(
            summary,
            "resume_decision\t{}\t{}\t{}",
            decision.gene, decision.action, decision.reason
        )?;
    }
    for gene in &execution.mergeable_genes {
        let gene_id = GeneId::parse(gene)?;
        let manifest_path = gene_dir_path(&args.output_root, &gene_id)?.join(MANIFEST_FILE_NAME);
        if !manifest_path.exists() {
            continue;
        }
        let manifest = read_run_manifest(&manifest_path)
            .with_context(|| format!("read input hashes for summary gene {gene:?}"))?;
        for input in &manifest.request.inputs {
            writeln!(
                summary,
                "input_sha256\t{}\t{}\t{}\t{}",
                gene, input.role, input.sha256, input.path
            )?;
        }
    }
    writeln!(
        summary,
        "elapsed_seconds\t{}",
        execution.elapsed.as_secs_f64()
    )?;
    summary.flush().context("flush batch summary")
}

/// Publish all report artifacts. Each file is committed by an atomic rename.
pub(super) fn publish_batch_report(
    context: BatchReportContext<'_>,
    execution: &mut ExecutionReport,
) -> anyhow::Result<PublishedBatchReport> {
    execution
        .resume_decisions
        .sort_by(|left, right| left.gene.cmp(&right.gene));
    execution.error_lines.sort();
    execution
        .downsample_records
        .sort_by(|left, right| left.gene.cmp(&right.gene));

    let summary_path = context
        .args
        .output_root
        .join(format!("{}_summary.txt", context.batch_file_prefix));
    let error_path = context
        .args
        .output_root
        .join(format!("{}_errors.txt", context.batch_file_prefix));
    let downsample_path = context
        .args
        .output_root
        .join(format!("{}_downsample.tsv", context.batch_file_prefix));
    for path in [&summary_path, &error_path, &downsample_path] {
        ensure_destination_within(&context.args.output_root, path)?;
    }

    atomic_write_with(&summary_path, |writer| {
        write_summary(writer, &context, execution)
    })?;
    let first_error = execution.error_lines.first().cloned();
    if execution.error_lines.is_empty() {
        if error_path.exists() {
            fs::remove_file(&error_path)
                .with_context(|| format!("remove stale error report {error_path:?}"))?;
        }
    } else {
        atomic_write_with(&error_path, |output| {
            for line in &execution.error_lines {
                writeln!(output, "{line}")?;
            }
            output.flush().context("flush batch error report")
        })?;
    }

    if execution.infrastructure_error_count() == 0 {
        for pair in execution.downsample_records.windows(2) {
            if pair[0].gene == pair[1].gene {
                anyhow::bail!(
                    "duplicate downsample state collected for gene {:?}",
                    pair[0].gene
                );
            }
        }
        if execution.downsample_records.is_empty() {
            if downsample_path.exists() {
                fs::remove_file(&downsample_path).with_context(|| {
                    format!("remove stale aggregate downsample state {downsample_path:?}")
                })?;
            }
        } else {
            atomic_write_with(&downsample_path, |writer| {
                write_downsample_records_to_writer(
                    writer,
                    &downsample_path,
                    &execution.downsample_records,
                )
            })
            .with_context(|| format!("write aggregate downsample state {downsample_path:?}"))?;
        }
    }

    Ok(PublishedBatchReport {
        summary_path,
        error_path,
        downsample_path,
        first_error,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normalized_options_include_runtime_and_scientific_sections() {
        let args = BatchRunOptions {
            cluster_mode: crate::flow::full::ClusterMode::Clusterj,
            prepare_reads: None,
            prepare_reference: None,
            prepare_prefix: None,
            prepare: crate::flow::config::PrepareConfig::default(),
            prepare_rejected_read_tracks: 0,
            input_root: PathBuf::from("input"),
            gene_list: None,
            output_root: PathBuf::from("output"),
            clustering: crate::flow::config::ClusteringConfig::default(),
            counting: crate::flow::config::CountingConfig::default(),
            runtime: crate::flow::config::RuntimeConfig {
                threads: 8,
                ..crate::flow::config::RuntimeConfig::default()
            },
            downsample: crate::flow::config::DownsampleConfig::default(),
        };
        let options = effective_options_json(&args, 3);
        assert_eq!(options["runtime"]["requested_threads"], 8);
        assert_eq!(options["runtime"]["effective_threads"], 3);
        assert!(options.get("clustering").is_some());
        assert!(options.get("counting").is_some());
        assert!(options.get("downsample").is_some());
    }
}
