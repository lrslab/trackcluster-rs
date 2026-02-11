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
    #[arg(long = "fraction-read", default_value_t = 0.01)]
    pub fraction_read: f64,

    /// Minimum fraction of reference span overlapping a read span
    #[arg(long = "fraction-ref", default_value_t = 0.05)]
    pub fraction_ref: f64,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let res = crate::flow::preparedir::prepare_dir_from_paths(
        &args.reads,
        &args.reference,
        &args.output_root,
        &args.prefix,
        AddGeneOpts {
            fraction_read: args.fraction_read,
            fraction_ref: args.fraction_ref,
        },
    )?;

    eprintln!(
        "preparedir: output_root={:?} genes={} dedup_reads={} novel_reads={}",
        args.output_root,
        res.genes.len(),
        res.dedup_reads,
        res.novel_reads
    );

    Ok(())
}
