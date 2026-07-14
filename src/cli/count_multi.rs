use std::path::{Path, PathBuf};

use anyhow::Context;

fn append_output_suffix(prefix: &Path, suffix: &str) -> PathBuf {
    let mut value = prefix.as_os_str().to_os_string();
    value.push(suffix);
    PathBuf::from(value)
}

#[derive(clap::Args, Debug)]
pub struct Args {
    /// Sample manifest TSV with columns: sample, reads, optional group
    #[arg(long = "manifest")]
    pub manifest: PathBuf,

    /// Reference BED (sorted recommended)
    #[arg(short = 'r', long = "reference")]
    pub reference: PathBuf,

    /// Isoform BED (from pooled cluster/flow)
    #[arg(short = 'i', long = "isoform")]
    pub isoform: PathBuf,

    /// Optional read-to-isoform TSV mapping (fast path; from pooled flow output)
    #[arg(long = "read-to-isoform")]
    pub read_to_isoform: Option<PathBuf>,

    /// How reads with multiple isoform candidates are counted: fractional or unique
    #[arg(long = "assignment-mode", default_value_t = crate::count::AssignmentMode::Unique)]
    pub assignment_mode: crate::count::AssignmentMode,

    /// Junction tolerance in bp used only for unique read-to-isoform assignment
    #[arg(
        long = "unique-assignment-junction-offset",
        default_value_t = crate::count::DEFAULT_UNIQUE_ASSIGNMENT_JUNCTION_OFFSET
    )]
    pub unique_assignment_junction_offset: u32,

    /// Output prefix for aggregate count, long, matrix, optional group, and unique-provenance files
    #[arg(short = 'o', long = "out")]
    pub out_prefix: PathBuf,
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

pub fn run(args: Args) -> anyhow::Result<()> {
    let unique_assignment_options = crate::count::UniqueAssignmentOptions {
        junction_offset: args.unique_assignment_junction_offset,
    };
    let sample_rows = crate::io::manifest::read_manifest_tsv(&args.manifest)?;
    let mapping_path = args
        .read_to_isoform
        .clone()
        .or_else(|| guess_mapping_path(&args.isoform));

    let count_csv = append_output_suffix(&args.out_prefix, ".isoform_count.csv");
    let long_tsv = append_output_suffix(&args.out_prefix, ".isoform_usage.long.tsv");
    let matrix_tsv = append_output_suffix(&args.out_prefix, ".isoform_counts.matrix.tsv");
    let group_tsv = append_output_suffix(&args.out_prefix, ".isoform_usage.group.tsv");
    let provenance = append_output_suffix(&args.out_prefix, ".unique_assignment.provenance.tsv");
    let mut inputs = vec![
        ("sample manifest input", args.manifest.as_path()),
        ("reference input", args.reference.as_path()),
        ("isoform input", args.isoform.as_path()),
    ];
    if let Some(path) = mapping_path.as_deref() {
        inputs.push(("read-to-isoform input", path));
    }
    for row in &sample_rows {
        inputs.push(("sample reads input", row.reads.as_path()));
    }
    super::ensure_distinct_inputs_and_outputs(
        &inputs,
        &[
            ("aggregate count output", count_csv.as_path()),
            ("long-format usage output", long_tsv.as_path()),
            ("count-matrix output", matrix_tsv.as_path()),
            ("group-usage output", group_tsv.as_path()),
            ("assignment-provenance output", provenance.as_path()),
        ],
    )?;

    let isoforms: Vec<crate::model::Transcript> = crate::io::bed::read_bed12(&args.isoform)?
        .collect::<Result<Vec<_>, crate::io::bed::BedError>>(
    )?;

    let outputs = if let Some(mapping_path) = mapping_path.as_ref() {
        let pairs = crate::count::read_read_to_isoform_tsv(mapping_path)
            .with_context(|| format!("read mapping {mapping_path:?}"))?;
        if args.assignment_mode == crate::count::AssignmentMode::Unique {
            let reads = crate::count::multi::read_tagged_sample_reads(&sample_rows)?;
            crate::count::multi::run_count_multi_from_read_to_isoform_unique_with_options(
                &sample_rows,
                &isoforms,
                &reads,
                &pairs,
                &args.out_prefix,
                unique_assignment_options,
            )?
        } else {
            crate::count::multi::run_count_multi_from_read_to_isoform(
                &sample_rows,
                &isoforms,
                &pairs,
                &args.out_prefix,
            )?
        }
    } else {
        let has_subreads = crate::count::has_embedded_subreads(&isoforms)?;
        if !has_subreads {
            anyhow::bail!(
                "count-multi: no --read-to-isoform provided and no mapping file found next to {:?}; \
this isoform BED does not embed read IDs (likely from --name2-mode coverage|none). \
Provide --read-to-isoform or re-run clustering with --name2-mode full.",
                args.isoform
            );
        }

        let refs: Vec<crate::model::Transcript> = crate::io::bed::read_bed12(&args.reference)
            .with_context(|| format!("open reference {:?}", args.reference))?
            .collect::<Result<Vec<_>, crate::io::bed::BedError>>()
            .with_context(|| format!("parse reference {:?}", args.reference))?;

        if args.assignment_mode == crate::count::AssignmentMode::Unique {
            let reads = crate::count::multi::read_tagged_sample_reads(&sample_rows)?;
            let pairs = crate::count::read_to_isoform_from_subreads(&isoforms, &refs)?;
            crate::count::multi::run_count_multi_from_read_to_isoform_unique_with_options(
                &sample_rows,
                &isoforms,
                &reads,
                &pairs,
                &args.out_prefix,
                unique_assignment_options,
            )?
        } else {
            crate::count::multi::run_count_multi(&sample_rows, &isoforms, &refs, &args.out_prefix)?
        }
    };

    eprintln!(
        "count-multi: count={:?} long={:?} matrix={:?} group={:?} provenance={:?}",
        outputs.count_csv,
        outputs.long_tsv,
        outputs.matrix_tsv,
        outputs.group_tsv,
        outputs.unique_assignment_provenance_tsv
    );

    Ok(())
}
