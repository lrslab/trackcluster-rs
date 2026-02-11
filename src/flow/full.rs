use std::ffi::OsString;
use std::fs;
use std::io::{BufRead, BufReader, Read, Write};
use std::panic::{self, AssertUnwindSafe};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::time::Instant;

use anyhow::Context;

use crate::annotate::addgene::AddGeneOpts;
use crate::annotate::desc::{describe, DescOpts};
use crate::count::{count_by_subreads, write_counts_csv};
use crate::flow::preparedir::{prepare_dir_from_paths, PrepareDirResult};

#[derive(Clone, Debug)]
pub struct BatchRunOptions {
    pub prepare_reads: Option<PathBuf>,
    pub prepare_reference: Option<PathBuf>,
    pub prepare_prefix: Option<String>,
    pub prepare_fraction_read: f64,
    pub prepare_fraction_ref: f64,
    pub input_root: PathBuf,
    pub gene_list: Option<PathBuf>,
    pub output_root: PathBuf,
    pub threads: usize,
    pub sw_score: i64,
    pub batch_size: usize,
    pub batch_rounds: usize,
    pub force: bool,
    pub progress_every: usize,
}

#[derive(Clone, Debug)]
pub struct BatchRunResult {
    pub prepared: Option<PrepareDirResult>,
    pub total_genes: usize,
    pub processed: usize,
    pub skipped: usize,
    pub errors: usize,
    pub elapsed_seconds: f64,
    pub summary_path: PathBuf,
    pub error_path: PathBuf,
}

#[derive(Clone, Debug)]
pub struct FullFlowOptions {
    pub reads: PathBuf,
    pub reference: PathBuf,
    pub output_root: PathBuf,
    pub prefix: String,
    pub threads: usize,
    pub sw_score: i64,
    pub batch_size: usize,
    pub batch_rounds: usize,
    pub prepare_fraction_read: f64,
    pub prepare_fraction_ref: f64,
    pub force: bool,
    pub progress_every: usize,
}

#[derive(Clone, Debug)]
pub struct FullFlowResult {
    pub batch: BatchRunResult,
    pub isoform_bed: PathBuf,
    pub unused_bed: PathBuf,
    pub count_csv: PathBuf,
    pub desc_prefix: PathBuf,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum GeneOutcome {
    Processed,
    Skipped,
}

fn read_gene_list(path: &Path) -> anyhow::Result<Vec<String>> {
    let file = fs::File::open(path).with_context(|| format!("open gene list {path:?}"))?;
    let reader = BufReader::new(file);
    let mut genes: Vec<String> = Vec::new();
    for line in reader.lines() {
        let line = line.with_context(|| format!("read gene list {path:?}"))?;
        let gene = line.trim();
        if gene.is_empty() || gene.starts_with('#') {
            continue;
        }
        genes.push(gene.to_owned());
    }
    genes.sort();
    genes.dedup();
    Ok(genes)
}

fn discover_genes(root: &Path) -> anyhow::Result<Vec<String>> {
    let mut genes: Vec<String> = Vec::new();
    for entry in fs::read_dir(root).with_context(|| format!("read_dir {root:?}"))? {
        let entry = entry.with_context(|| format!("read_dir entry in {root:?}"))?;
        let ty = entry
            .file_type()
            .with_context(|| format!("stat {:?}", entry.path()))?;
        if !ty.is_dir() {
            continue;
        }
        genes.push(entry.file_name().to_string_lossy().to_string());
    }
    genes.sort();
    genes.dedup();
    Ok(genes)
}

fn process_gene(gene: &str, args: &BatchRunOptions) -> anyhow::Result<GeneOutcome> {
    let gene_dir = args.input_root.join(gene);
    let reads = gene_dir.join(format!("{gene}_nano.bed"));
    let reference = gene_dir.join(format!("{gene}_gff.bed"));
    if !reads.exists() || !reference.exists() {
        return Ok(GeneOutcome::Skipped);
    }

    let reads_len = fs::metadata(&reads)
        .with_context(|| format!("stat reads {reads:?}"))?
        .len();
    if reads_len == 0 {
        return Ok(GeneOutcome::Skipped);
    }

    let out_dir = args.output_root.join(gene);
    fs::create_dir_all(&out_dir).with_context(|| format!("create {out_dir:?}"))?;

    let out_isoforms = out_dir.join(format!("{gene}_simple_coveragej.bed"));
    let out_unused = out_dir.join(format!("{gene}_unused.bed"));
    let out_mapping = out_dir.join(format!("{gene}_read_to_isoform.tsv"));

    if !args.force && out_isoforms.exists() {
        return Ok(GeneOutcome::Skipped);
    }

    let reads: Vec<crate::model::Transcript> = crate::io::bed::read_bed12(&reads)
        .with_context(|| format!("open reads {reads:?}"))?
        .collect::<Result<Vec<_>, crate::io::bed::BedError>>()
        .with_context(|| format!("parse reads {reads:?}"))?;

    let refs: Vec<crate::model::Transcript> = crate::io::bed::read_bed12(&reference)
        .with_context(|| format!("open reference {reference:?}"))?
        .collect::<Result<Vec<_>, crate::io::bed::BedError>>()
        .with_context(|| format!("parse reference {reference:?}"))?;

    let result = crate::cluster::clusterj::clusterj(
        &reads,
        Some(&refs),
        1,
        args.sw_score,
        args.batch_size,
        args.batch_rounds,
    );

    crate::cluster::output::write_isoforms_bed(&out_isoforms, &result.isoforms)
        .with_context(|| format!("write {out_isoforms:?}"))?;
    crate::cluster::output::write_isoforms_bed(&out_unused, &result.unused)
        .with_context(|| format!("write {out_unused:?}"))?;
    crate::cluster::output::write_read_to_isoform_tsv(&out_mapping, &result.read_to_isoform)
        .with_context(|| format!("write {out_mapping:?}"))?;

    Ok(GeneOutcome::Processed)
}

pub fn run_clusterj_batch(args: BatchRunOptions) -> anyhow::Result<BatchRunResult> {
    fs::create_dir_all(&args.output_root)
        .with_context(|| format!("create {:?}", args.output_root))?;

    let prepared = if args.prepare_reads.is_some()
        || args.prepare_reference.is_some()
        || args.prepare_prefix.is_some()
    {
        let reads = args
            .prepare_reads
            .as_ref()
            .context("--prepare-reads is required when using prepare")?;
        let reference = args
            .prepare_reference
            .as_ref()
            .context("--prepare-reference is required when using prepare")?;
        let prefix = args
            .prepare_prefix
            .as_deref()
            .context("--prepare-prefix is required when using prepare")?;

        eprintln!(
            "prepare: reads={:?} reference={:?} input_root={:?} prefix={}",
            reads, reference, args.input_root, prefix
        );
        let res = prepare_dir_from_paths(
            reads,
            reference,
            &args.input_root,
            prefix,
            AddGeneOpts {
                fraction_read: args.prepare_fraction_read,
                fraction_ref: args.prepare_fraction_ref,
            },
        )?;
        eprintln!(
            "prepare: genes={}, dedup_reads={}, novel_reads={}",
            res.genes.len(),
            res.dedup_reads,
            res.novel_reads
        );
        Some(res)
    } else {
        None
    };

    let mut genes = match &args.gene_list {
        Some(list) => read_gene_list(list)?,
        None => discover_genes(&args.input_root)?,
    };
    genes.sort();
    genes.dedup();

    let total = genes.len();
    if total == 0 {
        anyhow::bail!("no genes found (input-root={:?})", args.input_root);
    }

    let started = Instant::now();
    eprintln!(
        "clusterj-batch: {} genes, {} worker threads",
        total,
        args.threads.max(1)
    );

    let processed = Arc::new(AtomicUsize::new(0));
    let skipped = Arc::new(AtomicUsize::new(0));
    let errors = Arc::new(AtomicUsize::new(0));
    let done = Arc::new(AtomicUsize::new(0));
    let error_lines: Arc<Mutex<Vec<String>>> = Arc::new(Mutex::new(Vec::new()));

    let queue = Arc::new(Mutex::new(genes));
    let worker_count = args.threads.max(1).min(total);

    let mut handles = Vec::with_capacity(worker_count);
    for _ in 0..worker_count {
        let queue = Arc::clone(&queue);
        let processed = Arc::clone(&processed);
        let skipped = Arc::clone(&skipped);
        let errors = Arc::clone(&errors);
        let done = Arc::clone(&done);
        let error_lines = Arc::clone(&error_lines);
        let args = args.clone();

        handles.push(std::thread::spawn(move || loop {
            let gene = {
                let mut guard = queue.lock().expect("work queue poisoned");
                guard.pop()
            };
            let Some(gene) = gene else {
                break;
            };

            let outcome = match panic::catch_unwind(AssertUnwindSafe(|| process_gene(&gene, &args)))
            {
                Ok(Ok(outcome)) => outcome,
                Ok(Err(err)) => {
                    errors.fetch_add(1, Ordering::Relaxed);
                    if let Ok(mut guard) = error_lines.lock() {
                        guard.push(format!("{gene}\t{err}"));
                    }
                    GeneOutcome::Skipped
                }
                Err(payload) => {
                    errors.fetch_add(1, Ordering::Relaxed);
                    let msg = payload
                        .downcast_ref::<&str>()
                        .copied()
                        .map(str::to_owned)
                        .or_else(|| payload.downcast_ref::<String>().cloned())
                        .unwrap_or_else(|| "<non-string panic>".to_owned());
                    if let Ok(mut guard) = error_lines.lock() {
                        guard.push(format!("{gene}\tpanic\t{msg}"));
                    }
                    GeneOutcome::Skipped
                }
            };

            match outcome {
                GeneOutcome::Processed => {
                    processed.fetch_add(1, Ordering::Relaxed);
                }
                GeneOutcome::Skipped => {
                    skipped.fetch_add(1, Ordering::Relaxed);
                }
            }

            let done_now = done.fetch_add(1, Ordering::Relaxed) + 1;
            let every = args.progress_every.max(1);
            if done_now.is_multiple_of(every) || done_now == total {
                eprintln!(
                    "progress: {done_now}/{total} (processed={}, skipped={}, errors={})",
                    processed.load(Ordering::Relaxed),
                    skipped.load(Ordering::Relaxed),
                    errors.load(Ordering::Relaxed),
                );
            }
        }));
    }

    for handle in handles {
        if handle.join().is_err() {
            errors.fetch_add(1, Ordering::Relaxed);
        }
    }

    let elapsed = started.elapsed();
    let processed = processed.load(Ordering::Relaxed);
    let skipped = skipped.load(Ordering::Relaxed);
    let errors = errors.load(Ordering::Relaxed);

    eprintln!(
        "done: processed={processed}, skipped={skipped}, errors={errors}, elapsed={:?}",
        elapsed
    );

    let summary_path = args.output_root.join("clusterj_batch_summary.txt");
    let mut summary =
        fs::File::create(&summary_path).with_context(|| format!("write {summary_path:?}"))?;
    writeln!(summary, "input_root\t{:?}", args.input_root)?;
    writeln!(summary, "gene_list\t{:?}", args.gene_list)?;
    writeln!(summary, "output_root\t{:?}", args.output_root)?;
    writeln!(summary, "threads\t{}", args.threads)?;
    writeln!(summary, "sw_score\t{}", args.sw_score)?;
    writeln!(summary, "batch_size\t{}", args.batch_size)?;
    writeln!(summary, "batch_rounds\t{}", args.batch_rounds)?;
    writeln!(summary, "force\t{}", args.force)?;
    writeln!(summary, "total_genes\t{}", total)?;
    writeln!(summary, "processed\t{}", processed)?;
    writeln!(summary, "skipped\t{}", skipped)?;
    writeln!(summary, "errors\t{}", errors)?;
    writeln!(summary, "elapsed_seconds\t{}", elapsed.as_secs_f64())?;

    let error_path = args.output_root.join("clusterj_batch_errors.txt");
    let error_lines = error_lines
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner());
    if !error_lines.is_empty() {
        let mut errors_file =
            fs::File::create(&error_path).with_context(|| format!("write {error_path:?}"))?;
        for line in error_lines.iter() {
            writeln!(errors_file, "{line}")?;
        }
    } else if error_path.exists() {
        let _ = fs::remove_file(&error_path);
    }

    Ok(BatchRunResult {
        prepared,
        total_genes: total,
        processed,
        skipped,
        errors,
        elapsed_seconds: elapsed.as_secs_f64(),
        summary_path,
        error_path,
    })
}

fn append_suffix(path: &Path, suffix: &str) -> PathBuf {
    let mut out: OsString = path.as_os_str().to_os_string();
    out.push(suffix);
    PathBuf::from(out)
}

fn merge_files(inputs: &[PathBuf], out: &Path) -> anyhow::Result<()> {
    let mut writer = std::io::BufWriter::new(
        fs::File::create(out).with_context(|| format!("create merged output {out:?}"))?,
    );
    for input in inputs {
        let mut bytes: Vec<u8> = Vec::new();
        let mut reader = std::io::BufReader::new(
            fs::File::open(input).with_context(|| format!("open {input:?}"))?,
        );
        reader
            .read_to_end(&mut bytes)
            .with_context(|| format!("read {input:?}"))?;
        if bytes.is_empty() {
            continue;
        }
        writer
            .write_all(&bytes)
            .with_context(|| format!("append {input:?} into {out:?}"))?;
        if !bytes.ends_with(b"\n") {
            writer
                .write_all(b"\n")
                .with_context(|| format!("final newline after {input:?}"))?;
        }
    }
    Ok(())
}

fn merge_gene_outputs(
    output_root: &Path,
    genes: &[String],
    per_gene_suffix: &str,
    merged_out: &Path,
) -> anyhow::Result<()> {
    let mut paths: Vec<PathBuf> = Vec::new();
    for gene in genes {
        let path = output_root
            .join(gene)
            .join(format!("{gene}{per_gene_suffix}"));
        if path.exists() {
            paths.push(path);
        }
    }
    paths.sort();
    merge_files(&paths, merged_out)?;
    Ok(())
}

fn run_count_and_desc(
    isoform_bed: &Path,
    reference_bed: &Path,
    count_csv: &Path,
    desc_prefix: &Path,
) -> anyhow::Result<()> {
    let isoforms: Vec<crate::model::Transcript> = crate::io::bed::read_bed12(isoform_bed)
        .with_context(|| format!("open isoform {isoform_bed:?}"))?
        .collect::<Result<Vec<_>, crate::io::bed::BedError>>()
        .with_context(|| format!("parse isoform {isoform_bed:?}"))?;
    let refs: Vec<crate::model::Transcript> = crate::io::bed::read_bed12(reference_bed)
        .with_context(|| format!("open reference {reference_bed:?}"))?
        .collect::<Result<Vec<_>, crate::io::bed::BedError>>()
        .with_context(|| format!("parse reference {reference_bed:?}"))?;

    let counts = count_by_subreads(&isoforms, &refs);
    write_counts_csv(count_csv, &counts)
        .with_context(|| format!("write count csv {count_csv:?}"))?;

    let desc = describe(&isoforms, &refs, DescOpts::default());

    let desc_path = append_suffix(desc_prefix, "_desc.txt");
    let mut writer = std::io::BufWriter::new(
        fs::File::create(&desc_path).with_context(|| format!("write {desc_path:?}"))?,
    );
    for row in &desc.desc_rows {
        writeln!(
            &mut writer,
            "{}\t{}\t{}\t{}\t{}",
            row.isoform_id, row.ref_id, row.gene, row.miss, row.extra
        )?;
    }

    let class4_path = append_suffix(desc_prefix, "_class4.txt");
    let mut writer = std::io::BufWriter::new(
        fs::File::create(&class4_path).with_context(|| format!("write {class4_path:?}"))?,
    );
    for row in &desc.class4_rows {
        writeln!(&mut writer, "{}\t{}", row.isoform_id, row.class)?;
    }

    let fusion_path = append_suffix(desc_prefix, "_fusion.txt");
    let mut writer = std::io::BufWriter::new(
        fs::File::create(&fusion_path).with_context(|| format!("write {fusion_path:?}"))?,
    );
    for row in &desc.fusion_rows {
        writeln!(&mut writer, "{}\t{}", row.isoform_id, row.genes.join(";"))?;
    }

    let class12_path = append_suffix(desc_prefix, "_class12.txt");
    let mut writer = std::io::BufWriter::new(
        fs::File::create(&class12_path).with_context(|| format!("write {class12_path:?}"))?,
    );
    for row in &desc.class12_rows {
        writeln!(&mut writer, "{}\t{}", row.isoform_id, row.class)?;
    }

    Ok(())
}

pub fn run_full_flow(opts: FullFlowOptions) -> anyhow::Result<FullFlowResult> {
    fs::create_dir_all(&opts.output_root)
        .with_context(|| format!("create {:?}", opts.output_root))?;

    let gene_list = opts.output_root.join(format!("{}_gene.txt", opts.prefix));
    let batch = run_clusterj_batch(BatchRunOptions {
        prepare_reads: Some(opts.reads.clone()),
        prepare_reference: Some(opts.reference.clone()),
        prepare_prefix: Some(opts.prefix.clone()),
        prepare_fraction_read: opts.prepare_fraction_read,
        prepare_fraction_ref: opts.prepare_fraction_ref,
        input_root: opts.output_root.clone(),
        gene_list: Some(gene_list.clone()),
        output_root: opts.output_root.clone(),
        threads: opts.threads,
        sw_score: opts.sw_score,
        batch_size: opts.batch_size,
        batch_rounds: opts.batch_rounds,
        force: opts.force,
        progress_every: opts.progress_every,
    })?;

    let genes = read_gene_list(&gene_list).with_context(|| format!("read {:?}", gene_list))?;
    let isoform_bed = opts
        .output_root
        .join(format!("{}_isoform.bed", opts.prefix));
    let unused_bed = opts.output_root.join(format!("{}_unused.bed", opts.prefix));
    let count_csv = opts
        .output_root
        .join(format!("{}_isoform_count.csv", opts.prefix));
    let desc_prefix = opts.output_root.join(&opts.prefix);

    eprintln!("flow: merge isoforms -> {:?}", isoform_bed);
    merge_gene_outputs(
        &opts.output_root,
        &genes,
        "_simple_coveragej.bed",
        &isoform_bed,
    )?;

    eprintln!("flow: merge unused -> {:?}", unused_bed);
    merge_gene_outputs(&opts.output_root, &genes, "_unused.bed", &unused_bed)?;

    eprintln!("flow: count + desc");
    run_count_and_desc(&isoform_bed, &opts.reference, &count_csv, &desc_prefix)?;

    Ok(FullFlowResult {
        batch,
        isoform_bed,
        unused_bed,
        count_csv,
        desc_prefix,
    })
}
