pub mod addgene;
pub mod bam2bigg;
pub mod cluster;
pub mod clusterj;
pub mod count;
pub mod count_multi;
pub mod desc;
pub mod export;
pub mod flow;
pub mod gff2bigg;
pub mod mod_aggregate;
pub mod mod_contrast;
pub mod mod_import_dorado;
pub mod mod_import_m6anet;
pub mod mod_site_summary;
pub mod mod_subsample;
pub mod preparedir;
pub mod validate_bed;

use clap::{Parser, Subcommand};

fn lexical_absolute(path: &std::path::Path) -> anyhow::Result<std::path::PathBuf> {
    let absolute = if path.is_absolute() {
        path.to_path_buf()
    } else {
        std::env::current_dir()
            .map_err(anyhow::Error::from)
            .map_err(|error| error.context("resolve current working directory"))?
            .join(path)
    };
    let mut normalized = std::path::PathBuf::new();
    for component in absolute.components() {
        match component {
            std::path::Component::Prefix(_) | std::path::Component::RootDir => {
                normalized.push(component.as_os_str());
            }
            std::path::Component::CurDir => {}
            std::path::Component::ParentDir => {
                normalized.pop();
            }
            std::path::Component::Normal(value) => normalized.push(value),
        }
    }
    Ok(normalized)
}

fn resolved_output_path(path: &std::path::Path) -> anyhow::Result<std::path::PathBuf> {
    use anyhow::Context as _;

    let absolute = lexical_absolute(path)?;
    let mut ancestor = absolute.as_path();
    loop {
        match std::fs::canonicalize(ancestor) {
            Ok(resolved_ancestor) => {
                let suffix = absolute
                    .strip_prefix(ancestor)
                    .expect("ancestor is derived from absolute path");
                return Ok(resolved_ancestor.join(suffix));
            }
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
                ancestor = ancestor.parent().with_context(|| {
                    format!("resolve nearest existing ancestor of output path {path:?}")
                })?;
            }
            Err(error) => {
                return Err(error).with_context(|| format!("resolve output path {path:?}"));
            }
        }
    }
}

fn existing_paths_refer_to_same_file(
    left: &std::path::Path,
    right: &std::path::Path,
) -> anyhow::Result<bool> {
    use anyhow::Context as _;

    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt as _;

        let metadata = |path: &std::path::Path| match std::fs::metadata(path) {
            Ok(metadata) => Ok(Some(metadata)),
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(None),
            Err(error) => Err(error).with_context(|| format!("stat path {path:?}")),
        };
        let Some(left_metadata) = metadata(left)? else {
            return Ok(false);
        };
        let Some(right_metadata) = metadata(right)? else {
            return Ok(false);
        };
        Ok(left_metadata.dev() == right_metadata.dev()
            && left_metadata.ino() == right_metadata.ino())
    }

    #[cfg(not(unix))]
    {
        let _ = (left, right);
        Ok(false)
    }
}

fn ensure_distinct_path_pair(
    left: &std::path::Path,
    left_label: &str,
    right: &std::path::Path,
    right_label: &str,
) -> anyhow::Result<()> {
    use anyhow::Context as _;

    let left_path = std::fs::canonicalize(left)
        .with_context(|| format!("resolve {left_label} path {left:?}"))?;
    let right_path = resolved_output_path(right)?;
    if left_path == right_path || existing_paths_refer_to_same_file(left, right)? {
        anyhow::bail!("{left_label} and {right_label} refer to the same file");
    }
    Ok(())
}

pub(crate) fn ensure_distinct_inputs_and_outputs(
    inputs: &[(&str, &std::path::Path)],
    outputs: &[(&str, &std::path::Path)],
) -> anyhow::Result<()> {
    for (input_label, input) in inputs {
        for (output_label, output) in outputs {
            ensure_distinct_path_pair(input, input_label, output, output_label)?;
        }
    }
    for (index, (left_label, left)) in outputs.iter().enumerate() {
        let left_path = resolved_output_path(left)?;
        for (right_label, right) in outputs.iter().skip(index + 1) {
            let right_path = resolved_output_path(right)?;
            if left_path == right_path || existing_paths_refer_to_same_file(left, right)? {
                anyhow::bail!("{left_label} and {right_label} refer to the same file");
            }
        }
    }
    Ok(())
}

pub(crate) fn ensure_distinct_input_output(
    input: &std::path::Path,
    output: &std::path::Path,
    input_kind: &str,
) -> anyhow::Result<()> {
    let input_label = format!("{input_kind} input");
    ensure_distinct_inputs_and_outputs(&[(input_label.as_str(), input)], &[("output", output)])
}

fn read_read_tracks(
    path: &std::path::Path,
    invalid_read_policy: crate::flow::config::InvalidReadPolicy,
) -> anyhow::Result<(
    Vec<crate::model::Transcript>,
    Vec<crate::io::bed::RejectedReadRecord>,
)> {
    use anyhow::Context as _;

    let mut reader =
        crate::io::bed::read_bed12(path).with_context(|| format!("open reads {path:?}"))?;
    let mut reads = Vec::new();
    loop {
        let next = match invalid_read_policy {
            crate::flow::config::InvalidReadPolicy::Skip => reader.next_recovering_read(),
            crate::flow::config::InvalidReadPolicy::Fail => reader.next_strict_read(),
        }
        .with_context(|| format!("parse reads {path:?}"))?;
        let Some(read) = next else {
            break;
        };
        reads.push(read);
    }
    Ok((reads, reader.take_rejected_reads()))
}

#[derive(Parser, Debug)]
#[command(
    name = "trackcluster",
    version,
    about = "Pure-Rust rewrite of TrackCluster (no Bedtools)"
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand, Debug)]
pub enum Commands {
    /// Validate BED12/bigGenePred structure, optionally reporting explicit legacy repairs.
    ValidateBed(validate_bed::Args),
    /// Cluster long reads against a reference catalog by splice-junction structure.
    Clusterj(clusterj::Args),
    /// Cluster long reads against a reference catalog by transcript overlap.
    Cluster(cluster::Args),
    /// Run preparation, per-gene clustering, merging, counting, and annotation.
    Flow(flow::Args),
    /// Count one sample against clustered isoforms.
    Count(count::Args),
    /// Count multiple samples from a manifest against clustered isoforms.
    #[command(name = "count-multi")]
    CountMulti(count_multi::Args),
    /// Assign gene annotations to reads by overlap with reference transcripts.
    #[command(name = "addgene")]
    AddGene(addgene::Args),
    /// Build per-gene input directories from reads and a reference catalog.
    #[command(name = "preparedir")]
    PrepareDir(preparedir::Args),
    /// Classify clustered isoforms relative to reference transcripts.
    Desc(desc::Args),
    /// Export a BED transcript catalog as GTF, GFF3, or SQANTI3 audit input.
    Export(export::Args),
    /// Convert genome-aligned BAM records to TrackCluster bigGenePred-compatible BED.
    #[command(name = "bam2bigg", visible_alias = "bam-to-bigg")]
    Bam2Bigg(bam2bigg::Args),
    /// Convert GFF3 or GTF transcript annotations to TrackCluster bigGenePred-compatible BED.
    #[command(name = "gff2bigg", visible_alias = "gff-to-bigg")]
    Gff2Bigg(gff2bigg::Args),
    /// Aggregate normalized read-site modification observations by unique isoform assignment.
    #[command(name = "mod-aggregate")]
    ModAggregate(mod_aggregate::Args),
    /// Summarize isoform/site audit rows into a deterministic per-site QC inventory.
    #[command(name = "mod-site-summary")]
    ModSiteSummary(mod_site_summary::Args),
    /// Import m6Anet RNA002 read probabilities into normalized genomic observations.
    #[command(name = "mod-import-m6anet")]
    ModImportM6anet(mod_import_m6anet::Args),
    /// Import genome-aligned Dorado/modBAM calls into normalized observations.
    #[command(name = "mod-import-dorado")]
    ModImportDorado(mod_import_dorado::Args),
    /// Calculate explicit shared-site modification effect sizes without inferential p-values.
    #[command(name = "mod-contrast")]
    ModContrast(mod_contrast::Args),
    /// Split one high-coverage sample into synchronized low-coverage pseudo-sample inputs.
    #[command(name = "mod-subsample")]
    ModSubsample(mod_subsample::Args),
}

pub fn run(cli: Cli) -> anyhow::Result<()> {
    match cli.command {
        Commands::ValidateBed(args) => validate_bed::run(args),
        Commands::Clusterj(args) => clusterj::run(args),
        Commands::Cluster(args) => cluster::run(args),
        Commands::Flow(args) => flow::run(args),
        Commands::Count(args) => count::run(args),
        Commands::CountMulti(args) => count_multi::run(args),
        Commands::AddGene(args) => addgene::run(args),
        Commands::PrepareDir(args) => preparedir::run(args),
        Commands::Desc(args) => desc::run(args),
        Commands::Export(args) => export::run(args),
        Commands::Bam2Bigg(args) => bam2bigg::run(args),
        Commands::Gff2Bigg(args) => gff2bigg::run(args),
        Commands::ModAggregate(args) => mod_aggregate::run(args),
        Commands::ModSiteSummary(args) => mod_site_summary::run(args),
        Commands::ModImportM6anet(args) => mod_import_m6anet::run(args),
        Commands::ModImportDorado(args) => mod_import_dorado::run(args),
        Commands::ModContrast(args) => mod_contrast::run(args),
        Commands::ModSubsample(args) => mod_subsample::run(args),
    }
}

pub fn run_from_env() -> anyhow::Result<()> {
    run(Cli::parse())
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::path::PathBuf;
    use std::time::{SystemTime, UNIX_EPOCH};

    use clap::Parser;

    use super::*;

    #[test]
    fn parses_standard_flags() {
        let cli = Cli::try_parse_from([
            "trackcluster",
            "clusterj",
            "-s",
            "reads.bed",
            "-r",
            "ref.bed",
            "-o",
            "out.bed",
            "-t",
            "4",
            "--sl-partial-5prime-offset",
            "16",
            "--sl-same-junction-5prime-offset",
            "26",
            "--sl-5prime-cluster-offset",
            "17",
            "--sl-5prime-min-support",
            "3",
            "--same-junction-3prime-offset",
            "60",
            "--3prime-cluster-offset",
            "12",
            "--3prime-min-support",
            "6",
            "--sw-score",
            "-1",
        ])
        .unwrap();

        match cli.command {
            Commands::Clusterj(args) => {
                assert_eq!(args.reads, PathBuf::from("reads.bed"));
                assert_eq!(args.reference, PathBuf::from("ref.bed"));
                assert_eq!(args.out, PathBuf::from("out.bed"));
                assert_eq!(args.threads, 4);
                assert_eq!(args.sl_partial_5prime_offset, Some(16));
                assert_eq!(args.sl_same_junction_5prime_offset, Some(26));
                assert_eq!(args.sl_5prime_cluster_offset, Some(17));
                assert_eq!(args.sl_5prime_min_support, Some(3));
                assert_eq!(args.same_junction_3prime_offset, Some(60));
                assert_eq!(args.three_prime_cluster_offset, Some(12));
                assert_eq!(args.three_prime_min_support, Some(6));
                assert_eq!(args.sw_score, -1);
            }
            _ => panic!("expected clusterj args"),
        }
    }

    #[test]
    fn clusterj_rna002_preset_expands_and_explicit_flags_override() {
        let cli = Cli::try_parse_from([
            "trackcluster",
            "clusterj",
            "-s",
            "reads.bed",
            "-r",
            "ref.bed",
            "--platform-preset",
            "rna002",
        ])
        .unwrap();

        match cli.command {
            Commands::Clusterj(args) => {
                let resolved = args.resolved_platform_options();
                assert_eq!(
                    args.platform_preset,
                    crate::cluster::clusterj::PlatformPreset::Rna002
                );
                assert_eq!(resolved.junction_correction.offset, 15);
                assert_eq!(resolved.junction_correction.min_support, 5);
                assert_eq!(resolved.sl_options.partial_five_prime_end_offset, 20);
                assert_eq!(resolved.sl_options.same_junction_five_prime_end_offset, 25);
                assert_eq!(resolved.sl_options.five_prime_cluster_offset, 20);
                assert_eq!(resolved.sl_options.min_five_prime_cluster_support, 2);
                assert_eq!(
                    resolved
                        .three_prime_options
                        .same_junction_three_prime_end_offset,
                    50
                );
                assert_eq!(resolved.three_prime_options.three_prime_cluster_offset, 15);
                assert_eq!(
                    resolved.three_prime_options.min_three_prime_cluster_support,
                    5
                );
            }
            _ => panic!("expected clusterj args"),
        }

        let cli = Cli::try_parse_from([
            "trackcluster",
            "clusterj",
            "-s",
            "reads.bed",
            "-r",
            "ref.bed",
            "--platform-preset",
            "rna002",
            "--junction-correction-offset",
            "12",
            "--sl-partial-5prime-offset",
            "18",
        ])
        .unwrap();

        match cli.command {
            Commands::Clusterj(args) => {
                let resolved = args.resolved_platform_options();
                assert_eq!(resolved.junction_correction.offset, 12);
                assert_eq!(resolved.sl_options.partial_five_prime_end_offset, 18);
                assert_eq!(resolved.sl_options.same_junction_five_prime_end_offset, 25);
                assert_eq!(resolved.sl_options.five_prime_cluster_offset, 20);
                assert_eq!(resolved.three_prime_options.three_prime_cluster_offset, 12);
            }
            _ => panic!("expected clusterj args"),
        }
    }

    #[test]
    fn clusterj_rna004_preset_uses_generic_cutoffs() {
        let cli = Cli::try_parse_from([
            "trackcluster",
            "clusterj",
            "-s",
            "reads.bed",
            "-r",
            "ref.bed",
            "--platform-preset",
            "rna004",
        ])
        .unwrap();

        match cli.command {
            Commands::Clusterj(args) => {
                let resolved = args.resolved_platform_options();
                assert_eq!(
                    args.platform_preset,
                    crate::cluster::clusterj::PlatformPreset::Rna004
                );
                assert_eq!(resolved.junction_correction.offset, 10);
                assert_eq!(resolved.junction_correction.min_support, 5);
                assert_eq!(resolved.sl_options.partial_five_prime_end_offset, 15);
                assert_eq!(resolved.sl_options.same_junction_five_prime_end_offset, 25);
                assert_eq!(resolved.sl_options.five_prime_cluster_offset, 15);
                assert_eq!(resolved.sl_options.min_five_prime_cluster_support, 2);
                assert_eq!(
                    resolved
                        .three_prime_options
                        .same_junction_three_prime_end_offset,
                    50
                );
                assert_eq!(resolved.three_prime_options.three_prime_cluster_offset, 10);
                assert_eq!(
                    resolved.three_prime_options.min_three_prime_cluster_support,
                    5
                );
            }
            _ => panic!("expected clusterj args"),
        }
    }

    #[test]
    fn parses_flow_flags() {
        let cli = Cli::try_parse_from([
            "trackcluster",
            "flow",
            "--cluster-mode",
            "cluster",
            "-s",
            "reads.bed",
            "-r",
            "ref.bed",
            "-o",
            "outdir",
            "--prefix",
            "sample",
            "--threads",
            "8",
            "--sw-score",
            "-1",
            "--unique-assignment-junction-offset",
            "7",
        ])
        .unwrap();

        match cli.command {
            Commands::Flow(args) => {
                assert_eq!(args.cluster_mode, crate::flow::full::ClusterMode::Cluster);
                assert_eq!(args.reads, Some(PathBuf::from("reads.bed")));
                assert_eq!(args.manifest, None);
                assert_eq!(args.reference, PathBuf::from("ref.bed"));
                assert_eq!(args.output_root, PathBuf::from("outdir"));
                assert_eq!(args.prefix, "sample");
                assert_eq!(args.threads, 8);
                assert_eq!(args.sw_score, -1);
                assert_eq!(args.unique_assignment_junction_offset, 7);
            }
            _ => panic!("expected flow args"),
        }
    }

    #[test]
    fn parses_cluster_flags() {
        let cli = Cli::try_parse_from([
            "trackcluster",
            "cluster",
            "-s",
            "reads.bed",
            "-r",
            "ref.bed",
            "--batch-size",
            "50",
            "--batch-rounds",
            "7",
            "--name2-mode",
            "none",
            "--sw-score",
            "-1",
        ])
        .unwrap();

        match cli.command {
            Commands::Cluster(args) => {
                assert_eq!(args.reads, PathBuf::from("reads.bed"));
                assert_eq!(args.reference, PathBuf::from("ref.bed"));
                assert_eq!(args.batch_size, 50);
                assert_eq!(args.batch_rounds, 7);
                assert_eq!(args.name2_mode, crate::cluster::clusterj::Name2Mode::None);
                assert_eq!(args.sw_score, -1);
            }
            _ => panic!("expected cluster args"),
        }
    }

    #[test]
    fn rejects_zero_batch_rounds_at_cli_boundary() {
        let error = Cli::try_parse_from([
            "trackcluster",
            "cluster",
            "-s",
            "reads.bed",
            "-r",
            "ref.bed",
            "--batch-rounds",
            "0",
        ])
        .expect_err("zero batch rounds must be rejected by clap");
        assert!(error
            .to_string()
            .contains("batch rounds must be at least 1"));
    }

    #[test]
    fn parses_flow_manifest_flags() {
        let cli = Cli::try_parse_from([
            "trackcluster",
            "flow",
            "--manifest",
            "samples.tsv",
            "-r",
            "ref.bed",
            "-o",
            "outdir",
            "--prefix",
            "pooled",
            "--emit-pooled-reads",
        ])
        .unwrap();

        match cli.command {
            Commands::Flow(args) => {
                assert_eq!(args.reads, None);
                assert_eq!(args.manifest, Some(PathBuf::from("samples.tsv")));
                assert_eq!(args.reference, PathBuf::from("ref.bed"));
                assert!(args.emit_pooled_reads);
            }
            _ => panic!("expected flow args"),
        }
    }

    #[test]
    fn parses_flow_count_only_without_reads_or_manifest() {
        let cli = Cli::try_parse_from([
            "trackcluster",
            "flow",
            "--count-only",
            "-r",
            "ref.bed",
            "-o",
            "outdir",
            "--prefix",
            "sample",
        ])
        .unwrap();

        match cli.command {
            Commands::Flow(args) => {
                assert!(args.count_only);
                assert_eq!(args.reads, None);
                assert_eq!(args.manifest, None);
                assert_eq!(args.reference, PathBuf::from("ref.bed"));
            }
            _ => panic!("expected flow args"),
        }
    }

    #[test]
    fn parses_strict_gene_error_policy_for_normal_flow() {
        let cli = Cli::try_parse_from([
            "trackcluster",
            "flow",
            "--reads",
            "reads.bed",
            "--reference",
            "ref.bed",
            "--output-root",
            "outdir",
            "--prefix",
            "sample",
            "--strict-gene-errors",
        ])
        .unwrap();

        match cli.command {
            Commands::Flow(args) => assert!(args.strict_gene_errors),
            _ => panic!("expected flow args"),
        }
    }

    #[test]
    fn parses_count_output_root_flags() {
        let cli = Cli::try_parse_from([
            "trackcluster",
            "count",
            "-r",
            "ref.bed",
            "-o",
            "outdir",
            "--prefix",
            "sample",
            "--cluster-mode",
            "cluster",
            "--unique-assignment-junction-offset",
            "9",
        ])
        .unwrap();

        match cli.command {
            Commands::Count(args) => {
                assert_eq!(args.reference, PathBuf::from("ref.bed"));
                assert_eq!(args.output_root, Some(PathBuf::from("outdir")));
                assert_eq!(args.prefix, Some("sample".to_owned()));
                assert_eq!(args.cluster_mode, crate::flow::full::ClusterMode::Cluster);
                assert_eq!(args.isoform, None);
                assert_eq!(args.unique_assignment_junction_offset, 9);
            }
            _ => panic!("expected count args"),
        }
    }

    #[test]
    fn parses_count_multi_flags() {
        let cli = Cli::try_parse_from([
            "trackcluster",
            "count-multi",
            "--manifest",
            "samples.tsv",
            "-r",
            "ref.bed",
            "-i",
            "isoform.bed",
            "-o",
            "out/prefix",
            "--unique-assignment-junction-offset",
            "11",
        ])
        .unwrap();

        match cli.command {
            Commands::CountMulti(args) => {
                assert_eq!(args.manifest, PathBuf::from("samples.tsv"));
                assert_eq!(args.reference, PathBuf::from("ref.bed"));
                assert_eq!(args.isoform, PathBuf::from("isoform.bed"));
                assert_eq!(args.out_prefix, PathBuf::from("out/prefix"));
                assert_eq!(args.assignment_mode, crate::count::AssignmentMode::Unique);
                assert_eq!(args.unique_assignment_junction_offset, 11);
            }
            _ => panic!("expected count-multi args"),
        }
    }

    fn distinct_path_test_dir(label: &str) -> PathBuf {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let root = std::env::temp_dir().join(format!(
            "trackcluster-distinct-path-{label}-{}-{nonce}",
            std::process::id()
        ));
        fs::create_dir_all(&root).unwrap();
        root
    }

    #[test]
    fn distinct_path_preflight_normalizes_nonexistent_parent_components() {
        let root = distinct_path_test_dir("lexical");
        let input = root.join("input.bed");
        fs::write(&input, "input\n").unwrap();
        let alias = root.join("new-directory").join("..").join("input.bed");
        let error = ensure_distinct_inputs_and_outputs(
            &[("input", input.as_path())],
            &[("output", alias.as_path())],
        )
        .unwrap_err();
        assert!(error.to_string().contains("refer to the same file"));
        fs::remove_dir_all(root).unwrap();
    }

    #[cfg(unix)]
    #[test]
    fn distinct_path_preflight_rejects_hard_linked_outputs() {
        let root = distinct_path_test_dir("hardlink");
        let input = root.join("input.bed");
        let left = root.join("left.tsv");
        let right = root.join("right.tsv");
        fs::write(&input, "input\n").unwrap();
        fs::write(&left, "previous\n").unwrap();
        fs::hard_link(&left, &right).unwrap();
        let error = ensure_distinct_inputs_and_outputs(
            &[("input", input.as_path())],
            &[
                ("left output", left.as_path()),
                ("right output", right.as_path()),
            ],
        )
        .unwrap_err();
        assert!(error.to_string().contains("refer to the same file"));
        fs::remove_dir_all(root).unwrap();
    }
}
