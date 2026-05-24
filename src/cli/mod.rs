pub mod addgene;
pub mod cluster;
pub mod clusterj;
pub mod count;
pub mod count_multi;
pub mod desc;
pub mod flow;
pub mod preparedir;
pub mod validate_bed;

use clap::{Parser, Subcommand};

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
    ValidateBed(validate_bed::Args),
    Clusterj(clusterj::Args),
    Cluster(cluster::Args),
    Flow(flow::Args),
    Count(count::Args),
    #[command(name = "count-multi")]
    CountMulti(count_multi::Args),
    #[command(name = "addgene")]
    AddGene(addgene::Args),
    #[command(name = "preparedir")]
    PrepareDir(preparedir::Args),
    Desc(desc::Args),
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
    }
}

pub fn run_from_env() -> anyhow::Result<()> {
    run(Cli::parse())
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

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
        ])
        .unwrap();

        match cli.command {
            Commands::Count(args) => {
                assert_eq!(args.reference, PathBuf::from("ref.bed"));
                assert_eq!(args.output_root, Some(PathBuf::from("outdir")));
                assert_eq!(args.prefix, Some("sample".to_owned()));
                assert_eq!(args.cluster_mode, crate::flow::full::ClusterMode::Cluster);
                assert_eq!(args.isoform, None);
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
        ])
        .unwrap();

        match cli.command {
            Commands::CountMulti(args) => {
                assert_eq!(args.manifest, PathBuf::from("samples.tsv"));
                assert_eq!(args.reference, PathBuf::from("ref.bed"));
                assert_eq!(args.isoform, PathBuf::from("isoform.bed"));
                assert_eq!(args.out_prefix, PathBuf::from("out/prefix"));
                assert_eq!(args.assignment_mode, crate::count::AssignmentMode::Unique);
            }
            _ => panic!("expected count-multi args"),
        }
    }
}
