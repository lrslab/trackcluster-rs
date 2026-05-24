use std::path::{Path, PathBuf};

use anyhow::Context;

#[derive(clap::Args, Debug)]
pub struct Args {
    /// Existing flow/cluster-batch output directory containing per-gene folders
    #[arg(short = 'o', long = "output-root")]
    pub output_root: Option<PathBuf>,

    /// Prefix for merged outputs like `<prefix>_isoform.bed`
    #[arg(long = "prefix")]
    pub prefix: Option<String>,

    /// Clustering algorithm used for the existing per-gene outputs: `clusterj` or `cluster`
    #[arg(long = "cluster-mode", default_value_t = crate::flow::full::ClusterMode::Clusterj)]
    pub cluster_mode: crate::flow::full::ClusterMode,

    /// Reads BED (legacy isoform-BED mode)
    #[arg(short = 's', long = "reads")]
    pub reads: Option<PathBuf>,

    /// Reference BED
    #[arg(short = 'r', long = "reference")]
    pub reference: PathBuf,

    /// Isoform BED (legacy low-level mode)
    #[arg(short = 'i', long = "isoform")]
    pub isoform: Option<PathBuf>,

    /// Optional read-to-isoform TSV mapping (legacy low-level mode)
    #[arg(long = "read-to-isoform")]
    pub read_to_isoform: Option<PathBuf>,

    /// How reads with multiple isoform candidates are counted: fractional or unique
    #[arg(long = "assignment-mode", default_value_t = crate::count::AssignmentMode::Unique)]
    pub assignment_mode: crate::count::AssignmentMode,

    /// Output CSV for legacy isoform-BED mode
    #[arg(long = "out")]
    pub out: Option<PathBuf>,
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

fn run_output_root_count(args: Args, output_root: PathBuf) -> anyhow::Result<()> {
    let has_legacy_inputs =
        args.reads.is_some() || args.isoform.is_some() || args.read_to_isoform.is_some();
    if has_legacy_inputs || args.out.is_some() {
        if args.prefix.is_none() {
            anyhow::bail!(
                "count: -o/--output-root is for recounting an output folder; \
use --out for legacy standalone isoform-BED count output"
            );
        }
        anyhow::bail!(
            "count: --output-root mode reads per-gene cluster outputs; \
do not combine it with --reads, --isoform, --read-to-isoform, or --out"
        );
    }
    let prefix = args
        .prefix
        .ok_or_else(|| anyhow::anyhow!("count: --output-root requires --prefix"))?;
    let assignment_mode = args.assignment_mode;

    let result = crate::flow::full::run_full_flow(crate::flow::full::FullFlowOptions {
        cluster_mode: args.cluster_mode,
        reads: None,
        manifest: None,
        reference: args.reference,
        output_root,
        prefix,
        threads: 1,
        sw_score: crate::cluster::clusterj::DEFAULT_SW_SCORE,
        batch_size: 500,
        batch_rounds: 100,
        name2_mode: crate::cluster::clusterj::Name2Mode::Coverage,
        platform_preset: crate::cluster::clusterj::PlatformPreset::Generic,
        junction_correction_options: crate::cluster::clusterj::JunctionCorrectionOptions::default(),
        sl_options: crate::cluster::clusterj::SlMergeOptions::default(),
        three_prime_options: crate::cluster::clusterj::ThreePrimeMergeOptions::default(),
        overlap_cutoff1: crate::cluster::cluster_overlap::DEFAULT_CUTOFF1,
        overlap_cutoff2: crate::cluster::cluster_overlap::DEFAULT_CUTOFF2,
        overlap_intron_weight: crate::cluster::cluster_overlap::DEFAULT_INTRON_WEIGHT,
        prepare_fraction_read: 0.01,
        prepare_fraction_ref: 0.05,
        assignment_mode,
        emit_pooled_reads: false,
        force: false,
        count_only: true,
        progress_every: 1000,
        heartbeat_seconds: 0,
        heartbeat_top: 5,
        downsample_genes: Vec::new(),
        max_reads_per_gene: 0,
        downsample_seed: 1,
    })?;

    if assignment_mode == crate::count::AssignmentMode::Unique {
        eprintln!(
            "count: isoform={:?} count={:?} unique_mapping={:?}",
            result.isoform_bed,
            result.count_csv,
            result.count_csv.with_file_name(format!(
                "{}_read_to_isoform.unique.tsv",
                result
                    .desc_prefix
                    .file_name()
                    .unwrap_or_default()
                    .to_string_lossy()
            ))
        );
    } else {
        eprintln!(
            "count: isoform={:?} count={:?}",
            result.isoform_bed, result.count_csv
        );
    }
    Ok(())
}

fn run_legacy_isoform_count(args: Args) -> anyhow::Result<()> {
    let reads_path = args.reads.as_ref().ok_or_else(|| {
        anyhow::anyhow!("count: provide --output-root/--prefix, or provide legacy --reads")
    })?;
    let isoform_path = args.isoform.as_ref().ok_or_else(|| {
        anyhow::anyhow!("count: provide --output-root/--prefix, or provide legacy --isoform")
    })?;
    if args.prefix.is_some() {
        anyhow::bail!("count: --prefix is only valid with --output-root");
    }

    let isoforms: Vec<crate::model::Transcript> = crate::io::bed::read_bed12(isoform_path)?
        .collect::<Result<Vec<_>, crate::io::bed::BedError>>(
    )?;

    let mapping_path = args
        .read_to_isoform
        .clone()
        .or_else(|| guess_mapping_path(isoform_path));

    let records = if let Some(mapping_path) = mapping_path.as_ref() {
        let mut pairs = crate::count::read_read_to_isoform_tsv(mapping_path)
            .with_context(|| format!("read mapping {mapping_path:?}"))?;
        if args.assignment_mode == crate::count::AssignmentMode::Unique {
            let reads: Vec<crate::model::Transcript> = crate::io::bed::read_bed12(reads_path)
                .with_context(|| format!("open reads {:?}", reads_path))?
                .collect::<Result<Vec<_>, crate::io::bed::BedError>>()
                .with_context(|| format!("parse reads {:?}", reads_path))?;
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
                isoform_path
            );
        }

        let refs: Vec<crate::model::Transcript> = crate::io::bed::read_bed12(&args.reference)?
            .collect::<Result<Vec<_>, crate::io::bed::BedError>>(
        )?;
        if args.assignment_mode == crate::count::AssignmentMode::Unique {
            let reads: Vec<crate::model::Transcript> = crate::io::bed::read_bed12(reads_path)
                .with_context(|| format!("open reads {:?}", reads_path))?
                .collect::<Result<Vec<_>, crate::io::bed::BedError>>()
                .with_context(|| format!("parse reads {:?}", reads_path))?;
            let pairs = crate::count::read_to_isoform_from_subreads(&isoforms, &refs);
            let pairs =
                crate::count::select_unique_best_read_to_isoform(&reads, &isoforms, &pairs)?;
            crate::count::count_by_read_to_isoform(&isoforms, &pairs)
        } else {
            crate::count::count_by_subreads(&isoforms, &refs)
        }
    };
    let out = args
        .out
        .unwrap_or_else(|| PathBuf::from("isoform_count.csv"));
    crate::count::write_counts_csv(&out, &records)?;

    Ok(())
}

pub fn run(args: Args) -> anyhow::Result<()> {
    if let Some(output_root) = args.output_root.clone() {
        return run_output_root_count(args, output_root);
    }
    run_legacy_isoform_count(args)
}
