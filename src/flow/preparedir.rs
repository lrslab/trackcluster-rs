use std::collections::{HashMap, HashSet};
use std::fs::{self, File, OpenOptions};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::Context;

use crate::annotate::addgene::{add_gene_no_dedup, dedup_longest_by_name, AddGeneOpts};
use crate::io::bed::{read_bed12, write_bed12_to_writer, BedError};
use crate::io::manifest::SampleRow;
use crate::model::{Strand, Transcript};
use crate::sample::{tagged_read_name, SAMPLE_DELIM};

const GENE_NAME_COL: usize = 5;
const BED_READ_BUFFER_BYTES: usize = 1024 * 1024;

// For large inputs, `preparedir` is implemented as a bucketed two-pass pipeline to bound memory.
// We keep the original in-memory implementation for smaller inputs to preserve deterministic
// per-gene sorting and the existing golden outputs.
const PREPARE_IN_MEMORY_MAX_BYTES: u64 = 64 * 1024 * 1024;

// Must be a power of two for `bucket = hash & (BUCKET_COUNT - 1)`.
const BUCKET_COUNT: usize = 256;

// Reference span index bin size.
const BIN_SIZE: u32 = 16_384;

// Max number of open per-gene read writers.
const GENE_WRITER_CACHE_CAPACITY: usize = 128;

#[derive(Clone, Debug)]
pub struct PrepareDirResult {
    pub genes: Vec<String>,
    pub dedup_reads: usize,
    pub novel_reads: usize,
}

fn gene_field(tx: &Transcript) -> &str {
    tx.extra_fields
        .get(GENE_NAME_COL)
        .map(|value| value.as_str())
        .unwrap_or(tx.name.as_str())
}

fn genes(tx: &Transcript) -> impl Iterator<Item = &str> {
    gene_field(tx)
        .split("||")
        .map(str::trim)
        .filter(|g| !g.is_empty() && *g != "none")
}

fn write_bed12_indices(
    path: &Path,
    transcripts: &[Transcript],
    indices: &[usize],
) -> Result<(), BedError> {
    let file = File::create(path).map_err(|source| BedError::IoWrite {
        path: path.to_path_buf(),
        source,
    })?;
    let mut writer = std::io::BufWriter::new(file);
    write_bed12_to_writer(&mut writer, indices.iter().map(|&idx| &transcripts[idx])).map_err(
        |source| BedError::IoWrite {
            path: path.to_path_buf(),
            source,
        },
    )?;
    Ok(())
}

fn sort_indices_by_coord(transcripts: &[Transcript], indices: &mut [usize]) {
    indices.sort_by(|&a, &b| {
        let left = &transcripts[a];
        let right = &transcripts[b];
        left.chrom
            .cmp(&right.chrom)
            .then_with(|| left.tx_start.cmp(&right.tx_start))
            .then_with(|| left.tx_end.cmp(&right.tx_end))
            .then_with(|| left.name.cmp(&right.name))
    });
}

fn prepare_dir(
    reads_raw: Vec<Transcript>,
    refs: Vec<Transcript>,
    output_root: &Path,
    prefix: &str,
    addgene_opts: AddGeneOpts,
) -> anyhow::Result<PrepareDirResult> {
    fs::create_dir_all(output_root).with_context(|| format!("create {output_root:?}"))?;

    // Step 1 (Python: list_to_dic + write prefix_dedup.bed)
    let reads_dedup = dedup_longest_by_name(&reads_raw);
    let dedup_path = output_root.join(format!("{prefix}_dedup.bed"));
    crate::io::bed::write_bed12(&dedup_path, reads_dedup.iter())
        .with_context(|| format!("write {dedup_path:?}"))?;

    // Step 2 (Python: intersect + tracklist_add_gene)
    let reads_annotated = add_gene_no_dedup(&reads_dedup, &refs, addgene_opts);

    // Group references by gene (Python: group_bigg_by_gene(bigg_ref))
    let mut ref_by_gene: HashMap<String, Vec<usize>> = HashMap::new();
    for (idx, tx) in refs.iter().enumerate() {
        for gene in genes(tx) {
            ref_by_gene.entry(gene.to_owned()).or_default().push(idx);
        }
    }

    // Group reads by gene (Python: group_bigg_by_gene(bigg_new))
    let mut reads_by_gene: HashMap<String, Vec<usize>> = HashMap::new();
    let mut novel_indices: Vec<usize> = Vec::new();
    for (idx, tx) in reads_annotated.iter().enumerate() {
        let field = gene_field(tx);
        if field.trim() == "none" || field.trim().is_empty() {
            novel_indices.push(idx);
            continue;
        }
        for gene in genes(tx) {
            reads_by_gene.entry(gene.to_owned()).or_default().push(idx);
        }
    }

    // Write novel reads file (Python: prefix_novel.bed)
    let novel_path = output_root.join(format!("{prefix}_novel.bed"));
    write_bed12_indices(&novel_path, &reads_annotated, &novel_indices)
        .with_context(|| format!("write {novel_path:?}"))?;

    // Materialize gene list (Python: name2file(genename_l, prefix_gene.txt))
    let mut genes: Vec<String> = reads_by_gene.keys().cloned().collect();
    genes.sort();

    let gene_list_path = output_root.join(format!("{prefix}_gene.txt"));
    let mut gene_list = BufWriter::new(
        File::create(&gene_list_path).with_context(|| format!("write {gene_list_path:?}"))?,
    );
    for gene in &genes {
        writeln!(gene_list, "{gene}")?;
    }

    // Create per-gene folders and write inputs (Python: write {gene}_gff.bed and {gene}_nano.bed)
    for gene in &genes {
        let gene_dir = output_root.join(gene);
        fs::create_dir_all(&gene_dir).with_context(|| format!("create {gene_dir:?}"))?;

        let mut ref_indices = ref_by_gene.get(gene).cloned().unwrap_or_default();
        sort_indices_by_coord(&refs, &mut ref_indices);
        let ref_path = gene_dir.join(format!("{gene}_gff.bed"));
        write_bed12_indices(&ref_path, &refs, &ref_indices)
            .with_context(|| format!("write {ref_path:?}"))?;

        let mut read_indices = reads_by_gene.get(gene).cloned().unwrap_or_default();
        sort_indices_by_coord(&reads_annotated, &mut read_indices);
        let reads_path = gene_dir.join(format!("{gene}_nano.bed"));
        write_bed12_indices(&reads_path, &reads_annotated, &read_indices)
            .with_context(|| format!("write {reads_path:?}"))?;
    }

    Ok(PrepareDirResult {
        genes,
        dedup_reads: reads_dedup.len(),
        novel_reads: novel_indices.len(),
    })
}

pub fn prepare_dir_from_paths(
    reads_bed: &Path,
    reference_bed: &Path,
    output_root: &Path,
    prefix: &str,
    addgene_opts: AddGeneOpts,
) -> anyhow::Result<PrepareDirResult> {
    let refs: Vec<Transcript> = read_bed12(reference_bed)
        .with_context(|| format!("open reference {reference_bed:?}"))?
        .collect::<Result<Vec<_>, BedError>>()
        .with_context(|| format!("parse reference {reference_bed:?}"))?;

    let reads_len = fs::metadata(reads_bed)
        .with_context(|| format!("stat reads {reads_bed:?}"))?
        .len();

    if reads_len <= PREPARE_IN_MEMORY_MAX_BYTES {
        let reads_raw: Vec<Transcript> = read_bed12(reads_bed)
            .with_context(|| format!("open reads {reads_bed:?}"))?
            .collect::<Result<Vec<_>, BedError>>()
            .with_context(|| format!("parse reads {reads_bed:?}"))?;
        return prepare_dir(reads_raw, refs, output_root, prefix, addgene_opts);
    }

    prepare_dir_from_paths_bucketed(reads_bed, &refs, output_root, prefix, addgene_opts)
}

pub fn prepare_dir_from_manifest_rows(
    sample_rows: &[SampleRow],
    reference_bed: &Path,
    output_root: &Path,
    prefix: &str,
    addgene_opts: AddGeneOpts,
    pooled_reads_out: Option<&Path>,
) -> anyhow::Result<PrepareDirResult> {
    fs::create_dir_all(output_root).with_context(|| format!("create {output_root:?}"))?;

    let refs: Vec<Transcript> = read_bed12(reference_bed)
        .with_context(|| format!("open reference {reference_bed:?}"))?
        .collect::<Result<Vec<_>, BedError>>()
        .with_context(|| format!("parse reference {reference_bed:?}"))?;

    let total_reads_len: u64 = sample_rows
        .iter()
        .map(|row| {
            fs::metadata(&row.reads)
                .with_context(|| format!("stat reads {:?}", row.reads))
                .map(|meta| meta.len())
        })
        .collect::<anyhow::Result<Vec<_>>>()?
        .into_iter()
        .sum();

    if total_reads_len <= PREPARE_IN_MEMORY_MAX_BYTES {
        return prepare_dir_from_manifest_rows_in_memory(
            sample_rows,
            &refs,
            output_root,
            prefix,
            addgene_opts,
            pooled_reads_out,
        );
    }

    prepare_dir_from_manifest_rows_bucketed(
        sample_rows,
        &refs,
        output_root,
        prefix,
        addgene_opts,
        pooled_reads_out,
    )
}

fn prepare_dir_from_manifest_rows_in_memory(
    sample_rows: &[SampleRow],
    refs: &[Transcript],
    output_root: &Path,
    prefix: &str,
    addgene_opts: AddGeneOpts,
    pooled_reads_out: Option<&Path>,
) -> anyhow::Result<PrepareDirResult> {
    let mut pooled_writer = if let Some(path) = pooled_reads_out {
        Some((
            path.to_path_buf(),
            BufWriter::new(
                File::create(path)
                    .with_context(|| format!("create pooled reads output {path:?}"))?,
            ),
        ))
    } else {
        None
    };

    let mut reads_raw: Vec<Transcript> = Vec::new();
    for row in sample_rows {
        let reader =
            read_bed12(&row.reads).with_context(|| format!("open reads {:?}", row.reads))?;
        for record in reader {
            let mut tx = record.with_context(|| format!("parse reads {:?}", row.reads))?;
            tx.name = tagged_read_name(&row.sample, &tx.name);
            if let Some((path, writer)) = pooled_writer.as_mut() {
                write_bed12_to_writer(writer, std::iter::once(&tx))
                    .with_context(|| format!("write pooled reads {path:?}"))?;
            }
            reads_raw.push(tx);
        }
    }

    if let Some((path, writer)) = pooled_writer.as_mut() {
        writer
            .flush()
            .with_context(|| format!("flush pooled reads output {path:?}"))?;
    }

    prepare_dir(reads_raw, refs.to_vec(), output_root, prefix, addgene_opts)
}

fn prepare_dir_from_paths_bucketed(
    reads_bed: &Path,
    refs: &[Transcript],
    output_root: &Path,
    prefix: &str,
    addgene_opts: AddGeneOpts,
) -> anyhow::Result<PrepareDirResult> {
    fs::create_dir_all(output_root).with_context(|| format!("create {output_root:?}"))?;

    let (bucket_dir, bucket_paths) = create_bucket_paths(output_root, prefix)?;
    partition_reads_file_into_buckets(reads_bed, &bucket_paths)
        .with_context(|| format!("bucket reads {reads_bed:?}"))?;
    let res = prepare_dir_from_buckets(&bucket_paths, refs, output_root, prefix, addgene_opts);
    let _ = fs::remove_dir_all(&bucket_dir);
    res
}

fn prepare_dir_from_manifest_rows_bucketed(
    sample_rows: &[SampleRow],
    refs: &[Transcript],
    output_root: &Path,
    prefix: &str,
    addgene_opts: AddGeneOpts,
    pooled_reads_out: Option<&Path>,
) -> anyhow::Result<PrepareDirResult> {
    fs::create_dir_all(output_root).with_context(|| format!("create {output_root:?}"))?;

    let (bucket_dir, bucket_paths) = create_bucket_paths(output_root, prefix)?;
    partition_manifest_into_buckets(sample_rows, &bucket_paths, pooled_reads_out)
        .with_context(|| "bucket reads from manifest".to_owned())?;
    let res = prepare_dir_from_buckets(&bucket_paths, refs, output_root, prefix, addgene_opts);
    let _ = fs::remove_dir_all(&bucket_dir);
    res
}

fn create_bucket_paths(
    output_root: &Path,
    prefix: &str,
) -> anyhow::Result<(PathBuf, Vec<PathBuf>)> {
    let dir = output_root.join(format!(
        ".trackcluster_preparedir_tmp_{}_{}",
        prefix,
        std::process::id()
    ));
    if dir.exists() {
        fs::remove_dir_all(&dir).with_context(|| format!("remove stale temp dir {dir:?}"))?;
    }
    fs::create_dir_all(&dir).with_context(|| format!("create temp dir {dir:?}"))?;

    let mut paths: Vec<PathBuf> = Vec::with_capacity(BUCKET_COUNT);
    for bucket in 0..BUCKET_COUNT {
        paths.push(dir.join(format!("reads_bucket_{bucket:03}.bed")));
    }
    Ok((dir, paths))
}

fn partition_reads_file_into_buckets(
    reads_bed: &Path,
    bucket_paths: &[PathBuf],
) -> anyhow::Result<()> {
    let mut writers: Vec<BufWriter<File>> = Vec::with_capacity(bucket_paths.len());
    for path in bucket_paths {
        let file = File::create(path).with_context(|| format!("create bucket {path:?}"))?;
        writers.push(BufWriter::new(file));
    }

    let file = File::open(reads_bed).with_context(|| format!("open reads {reads_bed:?}"))?;
    let mut reader = BufReader::with_capacity(BED_READ_BUFFER_BYTES, file);
    let mut line_buf = String::new();
    let mut line_no: usize = 0;

    loop {
        line_buf.clear();
        let read_len = reader
            .read_line(&mut line_buf)
            .with_context(|| format!("read reads {reads_bed:?}"))?;
        if read_len == 0 {
            break;
        }
        line_no += 1;

        let line = line_buf.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let name = bed_line_name(line).with_context(|| {
            format!("reads {reads_bed:?}:{line_no}: expected at least 4 columns")
        })?;
        let bucket = bucket_for_hash(fnv1a_hash(name.as_bytes()));
        writers[bucket]
            .write_all(line.as_bytes())
            .with_context(|| format!("write bucket {:?}", bucket_paths[bucket]))?;
        writers[bucket]
            .write_all(b"\n")
            .with_context(|| format!("write bucket {:?}", bucket_paths[bucket]))?;
    }

    for (bucket, writer) in writers.iter_mut().enumerate() {
        writer
            .flush()
            .with_context(|| format!("flush bucket {:?}", bucket_paths[bucket]))?;
    }

    Ok(())
}

fn partition_manifest_into_buckets(
    sample_rows: &[SampleRow],
    bucket_paths: &[PathBuf],
    pooled_reads_out: Option<&Path>,
) -> anyhow::Result<()> {
    let mut writers: Vec<BufWriter<File>> = Vec::with_capacity(bucket_paths.len());
    for path in bucket_paths {
        let file = File::create(path).with_context(|| format!("create bucket {path:?}"))?;
        writers.push(BufWriter::new(file));
    }

    let mut pooled_writer = if let Some(path) = pooled_reads_out {
        Some((
            path.to_path_buf(),
            BufWriter::new(
                File::create(path)
                    .with_context(|| format!("create pooled reads output {path:?}"))?,
            ),
        ))
    } else {
        None
    };

    for row in sample_rows {
        let file = File::open(&row.reads).with_context(|| format!("open reads {:?}", row.reads))?;
        let mut reader = BufReader::with_capacity(BED_READ_BUFFER_BYTES, file);
        let mut line_buf = String::new();
        let mut line_no: usize = 0;

        loop {
            line_buf.clear();
            let read_len = reader
                .read_line(&mut line_buf)
                .with_context(|| format!("read reads {:?}", row.reads))?;
            if read_len == 0 {
                break;
            }
            line_no += 1;

            let line = line_buf.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            if line.contains('\t') {
                let mut it = line.splitn(5, '\t');
                let c1 = it.next().unwrap_or_default();
                let c2 = it.next().with_context(|| {
                    format!(
                        "reads {:?}:{line_no}: expected at least 4 columns",
                        row.reads
                    )
                })?;
                let c3 = it.next().with_context(|| {
                    format!(
                        "reads {:?}:{line_no}: expected at least 4 columns",
                        row.reads
                    )
                })?;
                let name = it.next().with_context(|| {
                    format!(
                        "reads {:?}:{line_no}: expected at least 4 columns",
                        row.reads
                    )
                })?;
                let rest = it.next();

                let hash = fnv1a_hash_tagged(&row.sample, name);
                let bucket = bucket_for_hash(hash);

                write!(
                    writers[bucket],
                    "{c1}\t{c2}\t{c3}\t{}{}{}",
                    row.sample, SAMPLE_DELIM, name
                )
                .with_context(|| format!("write bucket {:?}", bucket_paths[bucket]))?;
                if let Some(rest) = rest {
                    write!(writers[bucket], "\t{rest}")
                        .with_context(|| format!("write bucket {:?}", bucket_paths[bucket]))?;
                }
                writeln!(writers[bucket])
                    .with_context(|| format!("write bucket {:?}", bucket_paths[bucket]))?;

                if let Some((path, writer)) = pooled_writer.as_mut() {
                    write!(
                        writer,
                        "{c1}\t{c2}\t{c3}\t{}{}{}",
                        row.sample, SAMPLE_DELIM, name
                    )
                    .with_context(|| format!("write pooled reads {path:?}"))?;
                    if let Some(rest) = rest {
                        write!(writer, "\t{rest}")
                            .with_context(|| format!("write pooled reads {path:?}"))?;
                    }
                    writeln!(writer).with_context(|| format!("write pooled reads {path:?}"))?;
                }
            } else {
                let fields: Vec<&str> = line.split_whitespace().collect();
                if fields.len() < 4 {
                    anyhow::bail!(
                        "reads {:?}:{line_no}: expected at least 4 columns",
                        row.reads
                    );
                }
                let name = fields[3];
                let hash = fnv1a_hash_tagged(&row.sample, name);
                let bucket = bucket_for_hash(hash);

                for (idx, field) in fields.iter().enumerate() {
                    if idx > 0 {
                        write!(writers[bucket], "\t")
                            .with_context(|| format!("write bucket {:?}", bucket_paths[bucket]))?;
                    }
                    if idx == 3 {
                        write!(writers[bucket], "{}{}{}", row.sample, SAMPLE_DELIM, name)
                            .with_context(|| format!("write bucket {:?}", bucket_paths[bucket]))?;
                    } else {
                        write!(writers[bucket], "{field}")
                            .with_context(|| format!("write bucket {:?}", bucket_paths[bucket]))?;
                    }
                }
                writeln!(writers[bucket])
                    .with_context(|| format!("write bucket {:?}", bucket_paths[bucket]))?;

                if let Some((path, writer)) = pooled_writer.as_mut() {
                    for (idx, field) in fields.iter().enumerate() {
                        if idx > 0 {
                            write!(writer, "\t")
                                .with_context(|| format!("write pooled reads {path:?}"))?;
                        }
                        if idx == 3 {
                            write!(writer, "{}{}{}", row.sample, SAMPLE_DELIM, name)
                                .with_context(|| format!("write pooled reads {path:?}"))?;
                        } else {
                            write!(writer, "{field}")
                                .with_context(|| format!("write pooled reads {path:?}"))?;
                        }
                    }
                    writeln!(writer).with_context(|| format!("write pooled reads {path:?}"))?;
                }
            }
        }
    }

    for (bucket, writer) in writers.iter_mut().enumerate() {
        writer
            .flush()
            .with_context(|| format!("flush bucket {:?}", bucket_paths[bucket]))?;
    }
    if let Some((path, writer)) = pooled_writer.as_mut() {
        writer
            .flush()
            .with_context(|| format!("flush pooled reads output {path:?}"))?;
    }

    Ok(())
}

fn prepare_dir_from_buckets(
    bucket_paths: &[PathBuf],
    refs: &[Transcript],
    output_root: &Path,
    prefix: &str,
    addgene_opts: AddGeneOpts,
) -> anyhow::Result<PrepareDirResult> {
    fs::create_dir_all(output_root).with_context(|| format!("create {output_root:?}"))?;

    let dedup_path = output_root.join(format!("{prefix}_dedup.bed"));
    let novel_path = output_root.join(format!("{prefix}_novel.bed"));
    let mut dedup_writer =
        BufWriter::new(File::create(&dedup_path).with_context(|| format!("write {dedup_path:?}"))?);
    let mut novel_writer =
        BufWriter::new(File::create(&novel_path).with_context(|| format!("write {novel_path:?}"))?);

    let mut ref_by_gene: HashMap<String, Vec<usize>> = HashMap::new();
    for (idx, tx) in refs.iter().enumerate() {
        for gene in genes(tx) {
            ref_by_gene.entry(gene.to_owned()).or_default().push(idx);
        }
    }

    let ref_index = RefOverlapIndex::build(refs);
    let mut overlap_scratch = RefOverlapScratch::new(refs.len());
    let mut gene_name_buf: Vec<&str> = Vec::new();

    let mut genes_seen: HashSet<String> = HashSet::new();
    let mut gene_writer_cache = GeneWriterCache::new(GENE_WRITER_CACHE_CAPACITY);

    let mut dedup_reads_total: usize = 0;
    let mut novel_reads_total: usize = 0;

    for bucket_path in bucket_paths {
        let bucket_len = fs::metadata(bucket_path)
            .with_context(|| format!("stat bucket {bucket_path:?}"))?
            .len();
        if bucket_len == 0 {
            continue;
        }

        let mut dedup_reads: Vec<Transcript> = Vec::new();
        let mut pos: HashMap<String, usize> = HashMap::new();

        let reader = read_bed12(bucket_path)
            .with_context(|| format!("open reads bucket {bucket_path:?}"))?;
        for record in reader {
            let tx = record.with_context(|| format!("parse reads bucket {bucket_path:?}"))?;
            match pos.get(&tx.name).copied() {
                None => {
                    dedup_reads.push(tx);
                    let idx = dedup_reads.len() - 1;
                    pos.insert(dedup_reads[idx].name.clone(), idx);
                }
                Some(existing_idx) => {
                    let existing_len = exon_len(&dedup_reads[existing_idx]);
                    let candidate_len = exon_len(&tx);
                    if candidate_len > existing_len {
                        dedup_reads[existing_idx] = tx;
                    }
                }
            }
        }

        dedup_reads_total += dedup_reads.len();
        write_bed12_to_writer(&mut dedup_writer, dedup_reads.iter())
            .with_context(|| format!("append {dedup_path:?}"))?;

        for tx in &mut dedup_reads {
            ref_index.collect_gene_names(
                tx,
                refs,
                addgene_opts,
                &mut overlap_scratch,
                &mut gene_name_buf,
            );
            if !gene_name_buf.is_empty() {
                let joined = gene_name_buf.join("||");
                set_extra(tx, GENE_NAME_COL, joined);
            }

            let field = gene_field(tx);
            if field.trim() == "none" || field.trim().is_empty() {
                novel_reads_total += 1;
                write_bed12_to_writer(&mut novel_writer, std::iter::once(&*tx))
                    .with_context(|| format!("append {novel_path:?}"))?;
                continue;
            }

            for gene in genes(tx) {
                let writer = ensure_gene_writer(
                    gene,
                    &mut genes_seen,
                    &mut gene_writer_cache,
                    &ref_by_gene,
                    refs,
                    output_root,
                )?;
                write_bed12_to_writer(writer, std::iter::once(&*tx)).with_context(|| {
                    format!("write reads for gene {gene:?} (output_root={output_root:?})")
                })?;
            }
        }
    }

    dedup_writer
        .flush()
        .with_context(|| format!("flush {dedup_path:?}"))?;
    novel_writer
        .flush()
        .with_context(|| format!("flush {novel_path:?}"))?;
    gene_writer_cache.flush_all()?;

    let mut genes: Vec<String> = genes_seen.into_iter().collect();
    genes.sort();
    let gene_list_path = output_root.join(format!("{prefix}_gene.txt"));
    let mut gene_list = BufWriter::new(
        File::create(&gene_list_path).with_context(|| format!("write {gene_list_path:?}"))?,
    );
    for gene in &genes {
        writeln!(gene_list, "{gene}")?;
    }
    gene_list
        .flush()
        .with_context(|| format!("flush {gene_list_path:?}"))?;

    Ok(PrepareDirResult {
        genes,
        dedup_reads: dedup_reads_total,
        novel_reads: novel_reads_total,
    })
}

fn ensure_gene_writer<'a>(
    gene: &str,
    genes_seen: &mut HashSet<String>,
    writer_cache: &'a mut GeneWriterCache,
    ref_by_gene: &HashMap<String, Vec<usize>>,
    refs: &[Transcript],
    output_root: &Path,
) -> anyhow::Result<&'a mut BufWriter<File>> {
    let is_new = if genes_seen.contains(gene) {
        false
    } else {
        genes_seen.insert(gene.to_owned());
        true
    };

    let gene_dir = output_root.join(gene);
    if is_new {
        fs::create_dir_all(&gene_dir).with_context(|| format!("create {gene_dir:?}"))?;

        let mut ref_indices = ref_by_gene.get(gene).cloned().unwrap_or_default();
        sort_indices_by_coord(refs, &mut ref_indices);
        let ref_path = gene_dir.join(format!("{gene}_gff.bed"));
        write_bed12_indices(&ref_path, refs, &ref_indices)
            .with_context(|| format!("write {ref_path:?}"))?;
    }

    let reads_path = gene_dir.join(format!("{gene}_nano.bed"));
    writer_cache.get_or_open(gene, &reads_path, is_new)
}

fn exon_len(tx: &Transcript) -> u32 {
    tx.exons.iter().map(|exon| exon.len()).sum()
}

fn span_len(tx: &Transcript) -> u32 {
    tx.tx_end.get().saturating_sub(tx.tx_start.get())
}

fn span_overlap_len(a: &Transcript, b: &Transcript) -> u32 {
    let start = a.tx_start.get().max(b.tx_start.get());
    let end = a.tx_end.get().min(b.tx_end.get());
    end.saturating_sub(start)
}

fn set_extra(tx: &mut Transcript, idx: usize, value: String) {
    if tx.extra_fields.len() <= idx {
        tx.extra_fields.resize(idx + 1, "none".to_owned());
    }
    tx.extra_fields[idx] = value;
}

const FNV_OFFSET_BASIS: u64 = 14695981039346656037;
const FNV_PRIME: u64 = 1099511628211;

fn fnv1a_hash(bytes: &[u8]) -> u64 {
    let mut hash: u64 = FNV_OFFSET_BASIS;
    for &b in bytes {
        hash ^= b as u64;
        hash = hash.wrapping_mul(FNV_PRIME);
    }
    hash
}

fn fnv1a_hash_tagged(sample: &str, read_name: &str) -> u64 {
    let mut hash: u64 = FNV_OFFSET_BASIS;
    for &b in sample.as_bytes() {
        hash ^= b as u64;
        hash = hash.wrapping_mul(FNV_PRIME);
    }
    for &b in SAMPLE_DELIM.as_bytes() {
        hash ^= b as u64;
        hash = hash.wrapping_mul(FNV_PRIME);
    }
    for &b in read_name.as_bytes() {
        hash ^= b as u64;
        hash = hash.wrapping_mul(FNV_PRIME);
    }
    hash
}

fn bucket_for_hash(hash: u64) -> usize {
    debug_assert!(BUCKET_COUNT.is_power_of_two());
    (hash as usize) & (BUCKET_COUNT - 1)
}

fn bed_line_name(line: &str) -> Option<&str> {
    if line.contains('\t') {
        line.split('\t').nth(3)
    } else {
        line.split_whitespace().nth(3)
    }
}

struct BinnedPartitionIndex {
    max_bin: u32,
    bins: Vec<Vec<usize>>,
}

impl BinnedPartitionIndex {
    fn build(refs: &[Transcript], indices: &[usize]) -> Self {
        let max_end = indices
            .iter()
            .map(|&idx| refs[idx].tx_end.get())
            .max()
            .unwrap_or(0);
        let max_bin = max_end / BIN_SIZE;

        let mut bins: Vec<Vec<usize>> = vec![Vec::new(); (max_bin + 1) as usize];
        for &idx in indices {
            let start = refs[idx].tx_start.get();
            let end = refs[idx].tx_end.get();
            if start >= end {
                continue;
            }

            let start_bin = start / BIN_SIZE;
            let end_bin = end.saturating_sub(1) / BIN_SIZE;
            for bin in start_bin..=end_bin {
                bins[bin as usize].push(idx);
            }
        }

        Self { max_bin, bins }
    }

    fn collect_candidates(
        &self,
        query_start: u32,
        query_end: u32,
        seen: &mut [u32],
        stamp: u32,
        out: &mut Vec<usize>,
    ) {
        if query_start >= query_end {
            return;
        }

        let start_bin = query_start / BIN_SIZE;
        if start_bin > self.max_bin {
            return;
        }

        let end_bin = (query_end.saturating_sub(1) / BIN_SIZE).min(self.max_bin);
        for bin in start_bin..=end_bin {
            for &idx in &self.bins[bin as usize] {
                if seen[idx] != stamp {
                    seen[idx] = stamp;
                    out.push(idx);
                }
            }
        }
    }
}

#[derive(Default)]
struct ChromBinsBuilder {
    plus: Vec<usize>,
    minus: Vec<usize>,
    unknown: Vec<usize>,
}

struct ChromBins {
    plus: Option<BinnedPartitionIndex>,
    minus: Option<BinnedPartitionIndex>,
    unknown: Option<BinnedPartitionIndex>,
}

impl ChromBins {
    fn for_strand(&self, strand: Strand) -> Option<&BinnedPartitionIndex> {
        match strand {
            Strand::Plus => self.plus.as_ref(),
            Strand::Minus => self.minus.as_ref(),
            Strand::Unknown => self.unknown.as_ref(),
        }
    }
}

struct RefOverlapIndex {
    by_chrom: HashMap<String, ChromBins>,
}

impl RefOverlapIndex {
    fn build(refs: &[Transcript]) -> Self {
        let mut builders: HashMap<String, ChromBinsBuilder> = HashMap::new();
        for (idx, tx) in refs.iter().enumerate() {
            let builder = builders.entry(tx.chrom.clone()).or_default();
            match tx.strand {
                Strand::Plus => builder.plus.push(idx),
                Strand::Minus => builder.minus.push(idx),
                Strand::Unknown => builder.unknown.push(idx),
            }
        }

        let mut by_chrom: HashMap<String, ChromBins> = HashMap::with_capacity(builders.len());
        for (chrom, builder) in builders {
            let plus = if builder.plus.is_empty() {
                None
            } else {
                Some(BinnedPartitionIndex::build(refs, &builder.plus))
            };
            let minus = if builder.minus.is_empty() {
                None
            } else {
                Some(BinnedPartitionIndex::build(refs, &builder.minus))
            };
            let unknown = if builder.unknown.is_empty() {
                None
            } else {
                Some(BinnedPartitionIndex::build(refs, &builder.unknown))
            };
            by_chrom.insert(
                chrom,
                ChromBins {
                    plus,
                    minus,
                    unknown,
                },
            );
        }

        Self { by_chrom }
    }

    fn collect_gene_names<'a>(
        &self,
        read: &Transcript,
        refs: &'a [Transcript],
        opts: AddGeneOpts,
        scratch: &mut RefOverlapScratch,
        out: &mut Vec<&'a str>,
    ) {
        out.clear();

        let Some(chrom_bins) = self.by_chrom.get(read.chrom.as_str()) else {
            return;
        };
        let Some(index) = chrom_bins.for_strand(read.strand) else {
            return;
        };

        let start = read.tx_start.get();
        let end = read.tx_end.get();
        if start >= end {
            return;
        }

        let stamp = scratch.next_stamp();
        scratch.candidates.clear();
        index.collect_candidates(
            start,
            end,
            &mut scratch.seen,
            stamp,
            &mut scratch.candidates,
        );

        let read_len = span_len(read);
        if read_len == 0 {
            return;
        }

        for &ref_idx in &scratch.candidates {
            let reference = &refs[ref_idx];
            let overlap = span_overlap_len(read, reference);
            if overlap == 0 {
                continue;
            }

            let ref_len = span_len(reference);
            if ref_len == 0 {
                continue;
            }

            let overlap_f = overlap as f64;
            if overlap_f / (read_len as f64) < opts.fraction_read {
                continue;
            }
            if overlap_f / (ref_len as f64) < opts.fraction_ref {
                continue;
            }

            out.push(gene_field(reference));
        }

        out.sort_unstable();
        out.dedup();
    }
}

struct RefOverlapScratch {
    seen: Vec<u32>,
    stamp: u32,
    candidates: Vec<usize>,
}

impl RefOverlapScratch {
    fn new(ref_len: usize) -> Self {
        Self {
            seen: vec![0; ref_len],
            stamp: 0,
            candidates: Vec::new(),
        }
    }

    fn next_stamp(&mut self) -> u32 {
        self.stamp = self.stamp.wrapping_add(1);
        if self.stamp == 0 {
            self.seen.fill(0);
            self.stamp = 1;
        }
        self.stamp
    }
}

struct GeneWriterEntry {
    writer: BufWriter<File>,
    last_used: u64,
}

struct GeneWriterCache {
    capacity: usize,
    tick: u64,
    entries: HashMap<String, GeneWriterEntry>,
}

impl GeneWriterCache {
    fn new(capacity: usize) -> Self {
        Self {
            capacity: capacity.max(1),
            tick: 0,
            entries: HashMap::new(),
        }
    }

    fn get_or_open<'a>(
        &'a mut self,
        gene: &str,
        path: &Path,
        truncate: bool,
    ) -> anyhow::Result<&'a mut BufWriter<File>> {
        self.tick = self.tick.wrapping_add(1);
        if self.entries.contains_key(gene) {
            let tick = self.tick;
            let entry = self
                .entries
                .get_mut(gene)
                .expect("contains_key checked in branch");
            entry.last_used = tick;
            return Ok(&mut entry.writer);
        }

        if self.entries.len() >= self.capacity {
            let lru_key = self
                .entries
                .iter()
                .min_by_key(|(_, entry)| entry.last_used)
                .map(|(key, _)| key.clone());
            if let Some(key) = lru_key {
                if let Some(mut entry) = self.entries.remove(&key) {
                    entry
                        .writer
                        .flush()
                        .with_context(|| format!("flush evicted writer {key:?} -> {path:?}"))?;
                }
            }
        }

        let file = if truncate {
            OpenOptions::new()
                .create(true)
                .write(true)
                .truncate(true)
                .open(path)
        } else {
            OpenOptions::new().create(true).append(true).open(path)
        }
        .with_context(|| format!("open gene reads {gene:?} at {path:?}"))?;

        self.entries.insert(
            gene.to_owned(),
            GeneWriterEntry {
                writer: BufWriter::new(file),
                last_used: self.tick,
            },
        );

        Ok(&mut self
            .entries
            .get_mut(gene)
            .expect("inserted entry missing")
            .writer)
    }

    fn flush_all(&mut self) -> anyhow::Result<()> {
        for (gene, entry) in self.entries.iter_mut() {
            entry
                .writer
                .flush()
                .with_context(|| format!("flush gene writer {gene:?}"))?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::time::{SystemTime, UNIX_EPOCH};

    use super::*;

    fn fresh_temp_dir(prefix: &str) -> PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("time went backwards")
            .as_nanos();
        let mut dir = std::env::temp_dir();
        dir.push(format!(
            "trackcluster_rs_{}_{}_{}",
            prefix,
            std::process::id(),
            nanos
        ));
        fs::create_dir_all(&dir).expect("create temp dir");
        dir
    }

    #[test]
    fn preparedir_bucketed_creates_gene_folders_and_files() {
        let repo_root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        let reads = repo_root.join("tests/fixtures/reads.bed");
        let reference = repo_root.join("tests/fixtures/ref.bed");

        let out_root = fresh_temp_dir("preparedir_bucketed");
        let prefix = "sample";

        let refs: Vec<Transcript> = read_bed12(&reference)
            .unwrap()
            .collect::<Result<Vec<_>, BedError>>()
            .unwrap();

        let res = prepare_dir_from_paths_bucketed(
            &reads,
            &refs,
            &out_root,
            prefix,
            AddGeneOpts::default(),
        )
        .expect("preparedir bucketed run");

        assert_eq!(res.genes, vec!["GENEA".to_owned()]);

        let gene_list = out_root.join(format!("{prefix}_gene.txt"));
        assert!(gene_list.exists());
        assert_eq!(fs::read_to_string(gene_list).unwrap(), "GENEA\n");

        let dedup = out_root.join(format!("{prefix}_dedup.bed"));
        assert!(dedup.exists());
        let dedup_reads: Vec<Transcript> = read_bed12(&dedup)
            .unwrap()
            .collect::<Result<Vec<_>, BedError>>()
            .unwrap();
        assert_eq!(dedup_reads.len(), 1);
        assert_eq!(dedup_reads[0].name, "read_trunc");
        assert_eq!(
            dedup_reads[0].extra_fields.get(5).map(|s| s.as_str()),
            Some("none")
        );

        let novel = out_root.join(format!("{prefix}_novel.bed"));
        assert!(novel.exists());
        let novel_reads: Vec<Transcript> = read_bed12(&novel)
            .unwrap()
            .collect::<Result<Vec<_>, BedError>>()
            .unwrap();
        assert_eq!(novel_reads.len(), 0);

        let gene_dir = out_root.join("GENEA");
        let gene_reads_path = gene_dir.join("GENEA_nano.bed");
        let gene_ref_path = gene_dir.join("GENEA_gff.bed");
        assert!(gene_reads_path.exists());
        assert!(gene_ref_path.exists());

        let gene_reads: Vec<Transcript> = read_bed12(&gene_reads_path)
            .unwrap()
            .collect::<Result<Vec<_>, BedError>>()
            .unwrap();
        assert_eq!(gene_reads.len(), 1);
        assert_eq!(
            gene_reads[0].extra_fields.get(5).map(|s| s.as_str()),
            Some("GENEA")
        );

        let gene_refs: Vec<Transcript> = read_bed12(&gene_ref_path)
            .unwrap()
            .collect::<Result<Vec<_>, BedError>>()
            .unwrap();
        assert_eq!(gene_refs.len(), 2);
        for tx in &gene_refs {
            assert_eq!(tx.extra_fields.get(5).map(|s| s.as_str()), Some("GENEA"));
        }
    }

    #[test]
    fn preparedir_manifest_bucketed_writes_pooled_reads_with_sample_tags() {
        let repo_root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        let manifest = repo_root.join("tests/fixtures/samples.tsv");
        let reference = repo_root.join("tests/fixtures/ref.bed");

        let out_root = fresh_temp_dir("preparedir_manifest_bucketed");
        let prefix = "pooled";

        let sample_rows = crate::io::manifest::read_manifest_tsv(&manifest).unwrap();
        let refs: Vec<Transcript> = read_bed12(&reference)
            .unwrap()
            .collect::<Result<Vec<_>, BedError>>()
            .unwrap();
        let pooled_out = out_root.join(format!("{prefix}_pooled_reads.bed"));

        let res = prepare_dir_from_manifest_rows_bucketed(
            &sample_rows,
            &refs,
            &out_root,
            prefix,
            AddGeneOpts::default(),
            Some(&pooled_out),
        )
        .expect("preparedir manifest bucketed run");

        assert_eq!(res.genes, vec!["GENEA".to_owned()]);
        assert!(pooled_out.exists());
        let pooled = fs::read_to_string(pooled_out).expect("read pooled reads");
        assert!(pooled.contains("\tS1::read_s1\t"));
        assert!(pooled.contains("\tS2::read_s2\t"));
    }
}
