use std::path::{Path, PathBuf};

use anyhow::Context;

#[derive(clap::Args, Debug)]
#[command(
    group(clap::ArgGroup::new("polya_input").required(true).args(["bam", "nanopolish", "polya_manifest"])),
    group(clap::ArgGroup::new("single_polya_input").args(["bam", "nanopolish"]))
)]
pub struct Args {
    /// Dorado BAM carrying pt:i tags (single sample; exact, unprefixed read IDs).
    #[arg(long, requires = "sample")]
    pub bam: Option<PathBuf>,

    /// Nanopolish polya TSV (single sample; only qc_tag=PASS contributes lengths).
    #[arg(long, requires = "sample")]
    pub nanopolish: Option<PathBuf>,

    /// Sample label for --bam or --nanopolish.
    #[arg(long, requires = "single_polya_input")]
    pub sample: Option<String>,

    /// Optional experimental group for --bam or --nanopolish.
    #[arg(long, requires = "single_polya_input")]
    pub group: Option<String>,

    /// TSV with sample and bam or nanopolish, optional group; use sample::read IDs.
    #[arg(long = "polya-manifest")]
    pub polya_manifest: Option<PathBuf>,

    /// Final isoform BED/bigGenePred catalog.
    #[arg(long)]
    pub isoforms: PathBuf,

    /// Final unique read-to-isoform TSV (two columns, no header).
    #[arg(long = "read-to-isoform")]
    pub read_to_isoform: PathBuf,

    /// Output prefix for .isoform_polya.tsv, .read_polya.tsv, and .polya_qc.tsv.
    #[arg(short, long, visible_alias = "out-prefix")]
    pub out: PathBuf,
}

pub(crate) fn run_with_inputs(
    inputs: &[crate::polya::SampleInput],
    mode: crate::polya::ReadIdMode,
    isoforms: &Path,
    read_to_isoform: &Path,
    prefix: &Path,
    manifest: Option<&Path>,
) -> anyhow::Result<()> {
    let mut paths = vec![
        ("isoform input", isoforms),
        ("read-to-isoform input", read_to_isoform),
    ];
    if let Some(manifest) = manifest {
        paths.push(("poly(A) manifest input", manifest));
    }
    paths.extend(
        inputs
            .iter()
            .map(|input| ("poly(A) source input", input.source.path())),
    );
    let output_paths = crate::polya::output_paths(prefix);
    let outputs = output_paths
        .iter()
        .map(|path| ("poly(A) output", path.as_path()))
        .collect::<Vec<_>>();
    super::ensure_distinct_inputs_and_outputs(&paths, &outputs)?;
    let isoforms = crate::io::bed::read_bed12(isoforms)?
        .collect::<Result<Vec<_>, crate::io::bed::BedError>>()?;
    let pairs = crate::count::read_read_to_isoform_tsv(read_to_isoform)
        .with_context(|| format!("read unique mapping {read_to_isoform:?}"))?;
    crate::polya::aggregate(inputs, isoforms, pairs, mode)?.write(prefix)
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let (inputs, mode) = if let Some(manifest) = &args.polya_manifest {
        (
            crate::polya::read_manifest(manifest)?,
            crate::polya::ReadIdMode::SamplePrefixed,
        )
    } else {
        let sample = args
            .sample
            .context("--bam/--nanopolish requires --sample")?;
        crate::polya::validate_sample(&sample)?;
        (
            vec![crate::polya::SampleInput {
                sample,
                group: args.group,
                source: match (args.bam, args.nanopolish) {
                    (Some(path), None) => crate::polya::InputSource::DoradoBam(path),
                    (None, Some(path)) => crate::polya::InputSource::NanopolishTsv(path),
                    _ => anyhow::bail!("provide exactly one of --bam or --nanopolish"),
                },
            }],
            crate::polya::ReadIdMode::Raw,
        )
    };
    run_with_inputs(
        &inputs,
        mode,
        &args.isoforms,
        &args.read_to_isoform,
        &args.out,
        args.polya_manifest.as_deref(),
    )
}

pub(crate) fn prepare_flow(
    args: &super::flow::Args,
) -> anyhow::Result<Option<(Vec<crate::polya::SampleInput>, crate::polya::ReadIdMode)>> {
    crate::flow::path_key::SafePathComponent::parse("output prefix", &args.prefix)?;
    let polya_prefix = args.output_root.join(&args.prefix);
    let sample_rows = args
        .manifest
        .as_deref()
        .map(crate::io::manifest::read_manifest_tsv)
        .transpose()?;
    let polya_inputs = if let Some(manifest) = &args.polya_manifest {
        let mut inputs = crate::polya::read_manifest(manifest)?;
        crate::polya::match_flow_samples(
            &mut inputs,
            sample_rows
                .as_deref()
                .context("flow: --polya-manifest requires --manifest")?,
        )?;
        Some((inputs, crate::polya::ReadIdMode::SamplePrefixed))
    } else if let Some(source) = args
        .polya_bam
        .as_ref()
        .map(|path| crate::polya::InputSource::DoradoBam(path.clone()))
        .or_else(|| {
            args.polya_nanopolish
                .as_ref()
                .map(|path| crate::polya::InputSource::NanopolishTsv(path.clone()))
        })
    {
        let sample = args.polya_sample.clone().unwrap_or_else(|| {
            source
                .path()
                .file_stem()
                .map(|value| value.to_string_lossy().into_owned())
                .unwrap_or_else(|| "sample".to_owned())
        });
        crate::polya::validate_sample(&sample)?;
        Some((
            vec![crate::polya::SampleInput {
                sample,
                group: None,
                source,
            }],
            crate::polya::ReadIdMode::Raw,
        ))
    } else {
        None
    };
    if let Some((inputs, _)) = &polya_inputs {
        if args.assignment_mode != crate::count::AssignmentMode::Unique {
            anyhow::bail!("flow: poly(A) aggregation requires --assignment-mode unique");
        }
        let mut extra_inputs = inputs
            .iter()
            .map(|input| ("poly(A) source input", input.source.path()))
            .collect::<Vec<_>>();
        if let Some(manifest) = args.polya_manifest.as_deref() {
            extra_inputs.push(("poly(A) manifest input", manifest));
        }
        crate::flow::path_key::reject_external_inputs_in_output_root(
            &args.output_root,
            extra_inputs,
        )?;
    }
    let polya_outputs = crate::polya::output_paths(&polya_prefix);
    let polya_outputs = polya_outputs
        .iter()
        .map(|path| ("poly(A) output", path.as_path()))
        .collect::<Vec<_>>();
    let mut flow_inputs = vec![("reference input", args.reference.as_path())];
    if let Some(reads) = args.reads.as_deref() {
        flow_inputs.push(("reads input", reads));
    }
    if let Some(manifest) = args.manifest.as_deref() {
        flow_inputs.push(("sample manifest input", manifest));
    }
    if let Some(manifest) = args.mod_manifest.as_deref() {
        flow_inputs.push(("modification manifest input", manifest));
    }
    if let Some(fasta) = args.mod_reference_fasta.as_deref() {
        flow_inputs.push(("modification reference input", fasta));
    }
    if let Some(contrasts) = args.mod_contrasts.as_deref() {
        flow_inputs.push(("modification contrasts input", contrasts));
    }
    if let Some(samples) = &sample_rows {
        flow_inputs.extend(
            samples
                .iter()
                .map(|row| ("sample reads input", row.reads.as_path())),
        );
    }
    let mod_rows = args
        .mod_manifest
        .as_deref()
        .map(crate::io::mod_manifest::read_mod_manifest_tsv)
        .transpose()?;
    if let Some(rows) = &mod_rows {
        for row in rows {
            flow_inputs.push((
                "modification observations input",
                row.observations.as_path(),
            ));
            flow_inputs.push(("modification metadata input", row.assay_metadata.as_path()));
            if let Some(bam) = row.coverage_bam.as_deref() {
                flow_inputs.push(("modification coverage BAM input", bam));
            }
        }
    }
    super::ensure_distinct_inputs_and_outputs(&flow_inputs, &polya_outputs)?;
    Ok(polya_inputs)
}
