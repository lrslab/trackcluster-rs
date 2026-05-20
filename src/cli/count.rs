use std::path::{Path, PathBuf};

use anyhow::Context;

#[derive(clap::Args, Debug)]
pub struct Args {
    /// Reads BED (sorted recommended)
    #[arg(short = 's', long = "reads")]
    pub reads: PathBuf,

    /// Reference BED (sorted recommended)
    #[arg(short = 'r', long = "reference")]
    pub reference: PathBuf,

    /// Isoform BED (from cluster/clusterj)
    #[arg(short = 'i', long = "isoform")]
    pub isoform: PathBuf,

    /// Optional read-to-isoform TSV mapping (fast path; from clusterj/flow outputs)
    #[arg(long = "read-to-isoform")]
    pub read_to_isoform: Option<PathBuf>,

    /// How reads with multiple isoform candidates are counted: fractional or unique
    #[arg(long = "assignment-mode", default_value_t = crate::count::AssignmentMode::Unique)]
    pub assignment_mode: crate::count::AssignmentMode,

    /// Output CSV
    #[arg(short = 'o', long = "out", default_value = "isoform_count.csv")]
    pub out: PathBuf,
}

fn guess_mapping_path(isoform: &Path) -> Option<PathBuf> {
    let candidate = isoform.with_extension("read_to_isoform.tsv");
    if candidate.exists() {
        return Some(candidate);
    }

    let file_name = isoform.file_name()?.to_string_lossy();
    let prefix = file_name.strip_suffix("_isoform.bed")?;
    let candidate = isoform.with_file_name(format!("{prefix}_read_to_isoform.tsv"));
    candidate.exists().then_some(candidate)
}

pub fn run(_args: Args) -> anyhow::Result<()> {
    let isoforms: Vec<crate::model::Transcript> = crate::io::bed::read_bed12(&_args.isoform)?
        .collect::<Result<Vec<_>, crate::io::bed::BedError>>(
    )?;

    let mapping_path = _args
        .read_to_isoform
        .clone()
        .or_else(|| guess_mapping_path(&_args.isoform));

    let records = if let Some(mapping_path) = mapping_path.as_ref() {
        let mut pairs = crate::count::read_read_to_isoform_tsv(mapping_path)
            .with_context(|| format!("read mapping {mapping_path:?}"))?;
        if _args.assignment_mode == crate::count::AssignmentMode::Unique {
            let reads: Vec<crate::model::Transcript> = crate::io::bed::read_bed12(&_args.reads)
                .with_context(|| format!("open reads {:?}", _args.reads))?
                .collect::<Result<Vec<_>, crate::io::bed::BedError>>()
                .with_context(|| format!("parse reads {:?}", _args.reads))?;
            pairs = crate::count::select_unique_best_read_to_isoform(&reads, &isoforms, &pairs)?;
        }
        crate::count::count_by_read_to_isoform(&isoforms, &pairs)
    } else {
        let has_subreads = isoforms
            .iter()
            .any(|tx| !crate::count::parse_subreads(tx).is_empty());
        if !has_subreads {
            anyhow::bail!(
                "count: no --read-to-isoform provided and no mapping file found next to {:?}; \
this isoform BED does not embed read IDs (likely from --name2-mode coverage|none). \
Provide --read-to-isoform or re-run clustering with --name2-mode full.",
                _args.isoform
            );
        }

        let refs: Vec<crate::model::Transcript> = crate::io::bed::read_bed12(&_args.reference)?
            .collect::<Result<Vec<_>, crate::io::bed::BedError>>(
        )?;
        if _args.assignment_mode == crate::count::AssignmentMode::Unique {
            let reads: Vec<crate::model::Transcript> = crate::io::bed::read_bed12(&_args.reads)
                .with_context(|| format!("open reads {:?}", _args.reads))?
                .collect::<Result<Vec<_>, crate::io::bed::BedError>>()
                .with_context(|| format!("parse reads {:?}", _args.reads))?;
            let pairs = crate::count::read_to_isoform_from_subreads(&isoforms, &refs);
            let pairs =
                crate::count::select_unique_best_read_to_isoform(&reads, &isoforms, &pairs)?;
            crate::count::count_by_read_to_isoform(&isoforms, &pairs)
        } else {
            crate::count::count_by_subreads(&isoforms, &refs)
        }
    };
    crate::count::write_counts_csv(&_args.out, &records)?;

    Ok(())
}
