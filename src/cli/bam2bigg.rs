use std::io::{BufWriter, Write};
use std::path::PathBuf;

use anyhow::Context;

#[derive(clap::Args, Debug)]
pub struct Args {
    /// Input BAM alignment file
    #[arg(short = 'b', long = "bamfile", visible_alias = "input")]
    pub bamfile: PathBuf,

    /// Output TrackCluster bigGenePred-compatible BED path
    #[arg(
        short = 'o',
        long = "out",
        visible_alias = "output",
        default_value = "bigg.bed"
    )]
    pub out: PathBuf,

    /// Minimum MAPQ retained
    #[arg(
        short = 's',
        long = "score",
        visible_alias = "min-mapq",
        default_value_t = 30
    )]
    pub score: u8,

    /// Sample/group label; defaults to the BAM file stem
    #[arg(short = 'g', long = "group")]
    pub group: Option<String>,

    /// Include alignments marked secondary (0x100)
    #[arg(long)]
    pub include_secondary: bool,

    /// Include alignments marked supplementary (0x800)
    #[arg(long)]
    pub include_supplementary: bool,
}

fn default_group(path: &std::path::Path) -> String {
    path.file_stem()
        .map(|value| value.to_string_lossy().into_owned())
        .filter(|value| !value.trim().is_empty())
        .unwrap_or_else(|| "none".to_owned())
}

pub fn run(args: Args) -> anyhow::Result<()> {
    super::ensure_distinct_input_output(&args.bamfile, &args.out, "BAM")?;
    let options = crate::io::bam::BamToBiggOptions {
        min_mapq: args.score,
        include_secondary: args.include_secondary,
        include_supplementary: args.include_supplementary,
        group: args.group.unwrap_or_else(|| default_group(&args.bamfile)),
    };
    let summary = crate::flow::artifact_manifest::atomic_write_with(&args.out, |temporary| {
        let mut writer = BufWriter::new(temporary);
        let summary =
            crate::io::bam::write_bam_to_bed_writer(&args.bamfile, &options, &mut writer)?;
        writer
            .flush()
            .with_context(|| format!("flush temporary BAM conversion output {:?}", args.out))?;
        Ok(summary)
    })?;
    eprintln!(
        "bam2bigg: total={} records={} skipped_unmapped={} skipped_secondary={} skipped_supplementary={} skipped_below_mapq={}",
        summary.total_records,
        summary.converted_records,
        summary.skipped_unmapped,
        summary.skipped_secondary,
        summary.skipped_supplementary,
        summary.skipped_below_mapq
    );
    Ok(())
}
