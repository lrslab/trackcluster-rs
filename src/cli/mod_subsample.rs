use std::path::PathBuf;

use clap::Args as ClapArgs;

use crate::modification::subsample::{create_subsample_bundle, SubsampleMode, SubsampleOptions};

fn parse_positive_usize(value: &str) -> Result<usize, String> {
    let value = value
        .parse::<usize>()
        .map_err(|_| format!("invalid positive integer {value:?}"))?;
    if value == 0 {
        return Err("value must be at least 1".to_owned());
    }
    Ok(value)
}

fn parse_replicates(value: &str) -> Result<usize, String> {
    let value = parse_positive_usize(value)?;
    if value > 1000 {
        return Err("replicates must not exceed 1000".to_owned());
    }
    Ok(value)
}

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// TrackCluster sample manifest containing the high-coverage source sample.
    #[arg(long)]
    pub manifest: PathBuf,
    /// Final unique read-to-isoform assignment TSV.
    #[arg(long = "read-to-isoform")]
    pub read_to_isoform: PathBuf,
    /// Modification manifest containing source observations, metadata, and coverage BAM.
    #[arg(long = "mod-manifest")]
    pub mod_manifest: PathBuf,
    /// High-coverage source sample identifier.
    #[arg(long = "source-sample")]
    pub source_sample: String,
    /// Prefix for generated sample IDs: <prefix>_001, <prefix>_002, and so on.
    #[arg(long = "sample-prefix", default_value = "subsample")]
    pub sample_prefix: String,
    /// Number of generated pseudo-samples.
    #[arg(long, default_value_t = 4, value_parser = parse_replicates)]
    pub replicates: usize,
    /// Source read molecules selected per pseudo-sample.
    #[arg(long = "reads-per-sample", value_parser = parse_positive_usize)]
    pub reads_per_sample: usize,
    /// Sampling relationship: disjoint or independent.
    #[arg(long, default_value = "disjoint")]
    pub mode: SubsampleMode,
    /// Deterministic base seed.
    #[arg(long, default_value_t = 1)]
    pub seed: u64,
    /// New output directory containing ready-to-run manifests and synchronized inputs.
    #[arg(long = "out-dir")]
    pub out_dir: PathBuf,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let result = create_subsample_bundle(&SubsampleOptions {
        manifest: args.manifest,
        read_to_isoform: args.read_to_isoform,
        mod_manifest: args.mod_manifest,
        source_sample: args.source_sample,
        sample_prefix: args.sample_prefix,
        replicates: args.replicates,
        reads_per_sample: args.reads_per_sample,
        mode: args.mode,
        seed: args.seed,
        out_dir: args.out_dir,
    })?;
    eprintln!(
        "mod-subsample: available_reads={} samples={} assays={} output={}",
        result.available_reads,
        result.samples.len(),
        result.assays,
        result.out_dir.display()
    );
    Ok(())
}
