use std::path::PathBuf;

use crate::annotate::addgene::AddGeneOpts;

#[derive(clap::Args, Debug)]
pub struct Args {
    /// Reads BED
    #[arg(short = 's', long = "reads")]
    pub reads: PathBuf,

    /// Reference BED
    #[arg(short = 'r', long = "reference")]
    pub reference: PathBuf,

    /// Output directory (created if missing)
    #[arg(short = 'o', long = "output-root")]
    pub output_root: PathBuf,

    /// Prefix used for outputs like `<prefix>_gene.txt` and `<prefix>_dedup.bed`
    #[arg(long = "prefix")]
    pub prefix: String,

    /// Minimum fraction of read span overlapping a reference span
    #[arg(long = "fraction-read", default_value_t = 0.01, allow_hyphen_values = true, value_parser = crate::config::parse_unit_fraction)]
    pub fraction_read: f64,

    /// Minimum fraction of reference span overlapping a read span
    #[arg(long = "fraction-ref", default_value_t = 0.05, allow_hyphen_values = true, value_parser = crate::config::parse_unit_fraction)]
    pub fraction_ref: f64,

    /// Handling of malformed read tracks: skip only that track, or fail preparation
    #[arg(long = "invalid-read-policy", default_value_t = crate::flow::config::InvalidReadPolicy::Skip)]
    pub invalid_read_policy: crate::flow::config::InvalidReadPolicy,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let res = crate::flow::preparedir::prepare_dir_from_paths_with_policy(
        &args.reads,
        &args.reference,
        &args.output_root,
        &args.prefix,
        AddGeneOpts {
            fraction_read: args.fraction_read,
            fraction_ref: args.fraction_ref,
        },
        args.invalid_read_policy,
    )?;

    eprintln!(
        "preparedir: output_root={:?} genes={} dedup_reads={} novel_reads={} rejected_read_tracks={} rejected_reads={:?}",
        args.output_root,
        res.genes.len(),
        res.dedup_reads,
        res.novel_reads,
        res.rejected_read_tracks,
        res.rejected_reads_path
    );

    Ok(())
}
