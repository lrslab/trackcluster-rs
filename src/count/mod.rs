use std::collections::{BTreeMap, HashMap, HashSet};
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::model::{Coord, Interval, Strand, Transcript};

pub mod multi;

pub const DEFAULT_UNIQUE_ASSIGNMENT_JUNCTION_OFFSET: u32 = 15;
const UNIQUE_CATALOG_BIN_SIZE: u32 = 16_384;

/// Options that control unique read-to-isoform assignment.
///
/// Callers that persist selected mappings should record these options alongside
/// the output so that the assignment can be reproduced.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct UniqueAssignmentOptions {
    /// Maximum absolute difference allowed at each end of a matched intron.
    pub junction_offset: u32,
}

impl Default for UniqueAssignmentOptions {
    fn default() -> Self {
        Self {
            junction_offset: DEFAULT_UNIQUE_ASSIGNMENT_JUNCTION_OFFSET,
        }
    }
}

#[derive(Debug, Default, Clone, PartialEq)]
pub struct CountRecord {
    pub gene: String,
    pub isoform_id: String,
    pub count: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum AssignmentMode {
    Fractional,
    #[default]
    Unique,
}

impl std::fmt::Display for AssignmentMode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Fractional => f.write_str("fractional"),
            Self::Unique => f.write_str("unique"),
        }
    }
}

impl std::str::FromStr for AssignmentMode {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.eq_ignore_ascii_case("fractional") {
            return Ok(Self::Fractional);
        }
        if s.eq_ignore_ascii_case("unique") {
            return Ok(Self::Unique);
        }
        Err(format!(
            "invalid assignment mode {s:?}; expected one of: fractional, unique"
        ))
    }
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct AssignmentScore {
    chrom_mismatch: bool,
    strand_mismatch: bool,
    unmatched_read_junctions: usize,
    extra_isoform_junctions: usize,
    terminal_delta: u64,
    exon_symmetric_difference: u64,
    exon_count_delta: usize,
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct CatalogKey {
    gene: Option<String>,
    chrom: String,
    strand: Strand,
}

impl CatalogKey {
    fn new(gene: Option<&str>, chrom: &str, strand: Strand) -> Self {
        Self {
            gene: gene.map(str::to_owned),
            chrom: chrom.to_owned(),
            strand,
        }
    }
}

#[derive(Clone, Debug, Default)]
struct SpanCandidateIndex {
    bins: HashMap<u32, Vec<usize>>,
}

impl SpanCandidateIndex {
    fn build(isoforms: &[Transcript], indices: &[usize]) -> Self {
        let mut bins: HashMap<u32, Vec<usize>> = HashMap::new();
        for &idx in indices {
            let isoform = &isoforms[idx];
            let start = isoform.tx_start.get();
            let end = isoform.tx_end.get();
            if start >= end {
                continue;
            }

            let start_bin = start / UNIQUE_CATALOG_BIN_SIZE;
            let end_bin = end.saturating_sub(1) / UNIQUE_CATALOG_BIN_SIZE;
            for bin in start_bin..=end_bin {
                bins.entry(bin).or_default().push(idx);
            }
        }

        Self { bins }
    }

    fn collect_overlapping(
        &self,
        read: &Transcript,
        isoforms: &[Transcript],
        seen: &mut [u32],
        stamp: u32,
        out: &mut Vec<usize>,
    ) {
        let start = read.tx_start.get();
        let end = read.tx_end.get();
        if start >= end {
            return;
        }

        let start_bin = start / UNIQUE_CATALOG_BIN_SIZE;
        let end_bin = end.saturating_sub(1) / UNIQUE_CATALOG_BIN_SIZE;
        for bin in start_bin..=end_bin {
            let Some(indices) = self.bins.get(&bin) else {
                continue;
            };
            for &idx in indices {
                if seen[idx] == stamp {
                    continue;
                }
                let isoform = &isoforms[idx];
                if isoform.tx_start < read.tx_end && read.tx_start < isoform.tx_end {
                    seen[idx] = stamp;
                    out.push(idx);
                }
            }
        }
    }
}

pub(crate) fn parse_subreads(
    tx: &Transcript,
) -> Result<Vec<String>, crate::identity::IdentityError> {
    let Some(name2) = tx.metadata().name2() else {
        return Ok(Vec::new());
    };
    crate::identity::decode_name2(name2)
}

pub(crate) fn has_embedded_subreads(
    isoforms: &[Transcript],
) -> Result<bool, crate::identity::IdentityError> {
    for isoform in isoforms {
        if !parse_subreads(isoform)?.is_empty() {
            return Ok(true);
        }
    }
    Ok(false)
}

pub fn count_by_subreads(
    isoforms: &[Transcript],
    _references: &[Transcript],
) -> anyhow::Result<Vec<CountRecord>> {
    let mut pairs = Vec::new();
    for isoform in isoforms {
        for name in parse_subreads(isoform)? {
            pairs.push((name, isoform.name.clone()));
        }
    }
    count_by_read_to_isoform(isoforms, &pairs)
}

pub fn read_to_isoform_from_subreads(
    isoforms: &[Transcript],
    _references: &[Transcript],
) -> anyhow::Result<Vec<(String, String)>> {
    crate::identity::validate_isoform_ids(isoforms)?;
    let mut pairs = Vec::new();
    for isoform in isoforms {
        for read_name in parse_subreads(isoform)? {
            pairs.push((read_name, isoform.name.clone()));
        }
    }
    Ok(pairs)
}

pub fn read_read_to_isoform_tsv<P: AsRef<Path>>(
    path: P,
) -> Result<Vec<(String, String)>, std::io::Error> {
    let file = std::fs::File::open(path)?;
    let reader = BufReader::new(file);

    let mut pairs: Vec<(String, String)> = Vec::new();
    for (line_no, line) in reader.lines().enumerate() {
        let line = line?;
        if line.is_empty() {
            continue;
        }

        let Some((read_id, isoform_id)) = line.split_once('\t') else {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("invalid read_to_isoform line {}: {:?}", line_no + 1, line),
            ));
        };

        if read_id.is_empty() || isoform_id.is_empty() || isoform_id.contains('\t') {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "read_to_isoform line {} must contain exactly two non-empty TSV fields",
                    line_no + 1
                ),
            ));
        }

        pairs.push((read_id.to_owned(), isoform_id.to_owned()));
    }

    Ok(pairs)
}

fn transcript_exon_len(tx: &Transcript) -> u64 {
    tx.exons.iter().map(|exon| u64::from(exon.len())).sum()
}

fn same_read_alignment(left: &Transcript, right: &Transcript) -> bool {
    left.chrom == right.chrom
        && left.strand == right.strand
        && left.tx_start == right.tx_start
        && left.tx_end == right.tx_end
        && left.exons == right.exons
        && gene_name(left) == gene_name(right)
}

fn read_structures_by_name(reads: &[Transcript]) -> anyhow::Result<HashMap<&str, &Transcript>> {
    let mut reads_by_name: HashMap<&str, (usize, &Transcript)> = HashMap::new();
    for (index, read) in reads.iter().enumerate() {
        let read_id = read.name.as_str();
        if read_id.is_empty() {
            anyhow::bail!("read id must not be empty at reads index {index}");
        }
        match reads_by_name.get(read_id).copied() {
            Some((first_index, existing)) if !same_read_alignment(existing, read) => {
                anyhow::bail!(
                    "read id {read_id:?} has conflicting alignments at reads indices {first_index} and {index}; unique assignment requires one unambiguous structure per molecule (identical duplicate rows are allowed)"
                );
            }
            Some(_) => {}
            None => {
                reads_by_name.insert(read_id, (index, read));
            }
        }
    }
    Ok(reads_by_name
        .into_iter()
        .map(|(read_id, (_, read))| (read_id, read))
        .collect())
}

fn transcript_exon_overlap(left: &Transcript, right: &Transcript) -> u64 {
    let mut overlap = 0u64;
    for left_exon in &left.exons {
        for right_exon in &right.exons {
            overlap += u64::from(left_exon.overlap_len(*right_exon));
        }
    }
    overlap
}

fn gene_name(tx: &Transcript) -> Option<&str> {
    tx.metadata()
        .gene_id()
        .map(str::trim)
        .filter(|gene| !gene.is_empty() && *gene != "none")
}

fn gene_tokens(field: &str) -> impl Iterator<Item = &str> {
    field
        .split("||")
        .map(str::trim)
        .filter(|gene| !gene.is_empty() && *gene != "none")
}

fn transcript_gene_tokens(tx: &Transcript) -> impl Iterator<Item = &str> {
    gene_name(tx).into_iter().flat_map(gene_tokens)
}

/// Gene metadata is compatible when either side is genuinely unannotated, or
/// when the two annotated gene sets intersect. Prepared reads may carry a
/// `GENEA||GENEB` field after overlap assignment, while a per-gene catalog
/// record normally carries only one of those identifiers.
fn gene_metadata_compatible(left: &Transcript, right: &Transcript) -> bool {
    let left_field = gene_name(left);
    let right_field = gene_name(right);
    let left_has_gene = left_field.is_some_and(|field| gene_tokens(field).next().is_some());
    let right_has_gene = right_field.is_some_and(|field| gene_tokens(field).next().is_some());

    if !left_has_gene || !right_has_gene {
        return true;
    }

    let left_field = left_field.expect("checked above");
    let right_field = right_field.expect("checked above");
    gene_tokens(left_field)
        .any(|left_gene| gene_tokens(right_field).any(|right_gene| left_gene == right_gene))
}

fn junctions_match(left: Interval, right: Interval, offset: u32) -> bool {
    left.start.get().abs_diff(right.start.get()) <= offset
        && left.end.get().abs_diff(right.end.get()) <= offset
}

fn is_three_prime_terminal_exon(isoform: &Transcript, exon_idx: usize) -> bool {
    match isoform.strand {
        Strand::Plus | Strand::Unknown => exon_idx + 1 == isoform.exons.len(),
        Strand::Minus => exon_idx == 0,
    }
}

fn interval_contains(container: Interval, contained: Interval) -> bool {
    container.start <= contained.start && contained.end <= container.end
}

fn lower_bound_introns_by_start(introns: &[Interval], start: Coord) -> usize {
    let mut left = 0usize;
    let mut right = introns.len();
    while left < right {
        let mid = left + (right - left) / 2;
        if introns[mid].start < start {
            left = mid + 1;
        } else {
            right = mid;
        }
    }
    left
}

fn retained_regions_by_isoform(isoforms: &[Transcript]) -> Vec<Vec<Interval>> {
    let mut out = vec![Vec::new(); isoforms.len()];
    let mut domains: HashMap<CatalogKey, Vec<usize>> = HashMap::new();
    for (isoform_idx, isoform) in isoforms.iter().enumerate() {
        domains
            .entry(CatalogKey::new(
                gene_name(isoform),
                &isoform.chrom,
                isoform.strand,
            ))
            .or_default()
            .push(isoform_idx);
    }

    // Missing gene metadata is an explicit domain of its own. It may share
    // introns only with other unannotated isoforms on the same chromosome and
    // strand, never with annotated genes or another chromosome/strand.
    for domain_indices in domains.values() {
        let mut domain_introns = Vec::new();
        for &isoform_idx in domain_indices {
            domain_introns.extend(isoforms[isoform_idx].introns());
        }
        domain_introns.sort_unstable();
        domain_introns.dedup();

        if domain_introns.is_empty() {
            continue;
        }

        for &isoform_idx in domain_indices {
            let isoform = &isoforms[isoform_idx];
            let retained_regions = &mut out[isoform_idx];
            for exon in &isoform.exons {
                let mut idx = lower_bound_introns_by_start(&domain_introns, exon.start);
                while idx < domain_introns.len() && domain_introns[idx].start <= exon.end {
                    let intron = domain_introns[idx];
                    if !intron.is_empty() && interval_contains(*exon, intron) {
                        retained_regions.push(intron);
                    }
                    idx += 1;
                }
            }

            retained_regions.sort_unstable();
            retained_regions.dedup();
        }
    }

    out
}

fn interval_spans_position(interval: Interval, position: Coord) -> bool {
    interval.start < position && position < interval.end
}

fn read_supports_retained_region(read: &Transcript, region: Interval) -> bool {
    read.exons.iter().any(|read_exon| {
        interval_spans_position(*read_exon, region.start)
            || interval_spans_position(*read_exon, region.end)
    })
}

fn retained_regions_supported_by_read(read: &Transcript, retained_regions: &[Interval]) -> bool {
    let read_span = read_span(read);
    let mut covered_regions = retained_regions
        .iter()
        .copied()
        .filter(|region| region.overlap_len(read_span) > 0)
        .peekable();

    // A partial read cannot establish retention when it never reaches a
    // retained region. When it reaches several regions, every reached region
    // must have a boundary spanned by a read exon; support for one retained
    // intron must not mask a splice through another.
    covered_regions.peek().is_some()
        && covered_regions.all(|region| read_supports_retained_region(read, region))
}

/// Return a maximum-cardinality, minimum-distance ordered intron alignment.
///
/// Every read intron and isoform intron can occur in at most one pair. Closely
/// spaced microfeatures therefore remain distinct: there is deliberately no
/// rule allowing one intron to satisfy multiple counterparts. Introns are in
/// ascending genomic order for both strands, which preserves their relative
/// order even though biological traversal is reversed on the minus strand.
/// Exact score ties retain the previously established prefix alignment, making
/// ambiguous repeated offsets deterministic.
pub(crate) fn ordered_one_to_one_intron_matches(
    read_introns: &[Interval],
    isoform_introns: &[Interval],
    offset: u32,
) -> Vec<(usize, usize)> {
    crate::matching::ordered_one_to_one_matches_by(
        read_introns.len(),
        isoform_introns.len(),
        |read_idx, isoform_idx| {
            let read_intron = read_introns[read_idx];
            let isoform_intron = isoform_introns[isoform_idx];
            junctions_match(read_intron, isoform_intron, offset).then(|| {
                u64::from(read_intron.start.get().abs_diff(isoform_intron.start.get()))
                    + u64::from(read_intron.end.get().abs_diff(isoform_intron.end.get()))
            })
        },
    )
}

fn matching_intron_indices(
    read_introns: &[Interval],
    isoform_introns: &[Interval],
    offset: u32,
) -> (HashSet<usize>, HashSet<usize>) {
    let matches = ordered_one_to_one_intron_matches(read_introns, isoform_introns, offset);
    let matched_read = matches.iter().map(|(read_idx, _)| *read_idx).collect();
    let matched_isoform = matches
        .iter()
        .map(|(_, isoform_idx)| *isoform_idx)
        .collect();
    (matched_read, matched_isoform)
}

fn unmatched_isoform_introns_overlapped_by_read_exons(
    read: &Transcript,
    isoform_introns: &[Interval],
    matched_isoform_introns: &HashSet<usize>,
) -> usize {
    isoform_introns
        .iter()
        .enumerate()
        .filter(|(idx, isoform_intron)| {
            !matched_isoform_introns.contains(idx)
                && read
                    .exons
                    .iter()
                    .any(|read_exon| read_exon.overlap_len(**isoform_intron) > 0)
        })
        .count()
}

fn junctionless_read_candidate(
    read: &Transcript,
    isoform: &Transcript,
    isoform_introns: &[Interval],
) -> bool {
    if read.exons.len() != 1 {
        return false;
    }

    let read_exon = read.exons[0];
    if isoform_introns
        .iter()
        .any(|isoform_intron| read_exon.overlap_len(*isoform_intron) > 0)
    {
        return false;
    }

    if isoform_introns.is_empty() {
        return true;
    }

    isoform.exons.iter().enumerate().any(|(idx, isoform_exon)| {
        is_three_prime_terminal_exon(isoform, idx) && interval_contains(*isoform_exon, read_exon)
    })
}

fn catalog_assignment_candidate(
    read: &Transcript,
    isoform: &Transcript,
    junction_offset: u32,
) -> bool {
    if read.chrom != isoform.chrom || read.strand != isoform.strand {
        return false;
    }
    if !gene_metadata_compatible(read, isoform) {
        return false;
    }
    if transcript_exon_overlap(read, isoform) == 0 {
        return false;
    }

    let read_introns = read.introns();
    let isoform_introns = isoform.introns();
    if read_introns.is_empty() {
        return junctionless_read_candidate(read, isoform, &isoform_introns);
    }
    if isoform_introns.is_empty() {
        return false;
    }

    let (matched_read_introns, matched_isoform_introns) =
        matching_intron_indices(&read_introns, &isoform_introns, junction_offset);

    matched_read_introns.len() == read_introns.len()
        && unmatched_isoform_introns_overlapped_by_read_exons(
            read,
            &isoform_introns,
            &matched_isoform_introns,
        ) == 0
}

fn covered_unmatched_isoform_junctions(
    read: &Transcript,
    isoform_introns: &[Interval],
    matched_isoform_introns: &HashSet<usize>,
) -> usize {
    isoform_introns
        .iter()
        .enumerate()
        .filter(|(idx, isoform_intron)| {
            if matched_isoform_introns.contains(idx) {
                return false;
            }

            let inside_read_span =
                isoform_intron.start >= read.tx_start && isoform_intron.end <= read.tx_end;
            let overlaps_read_exon = read
                .exons
                .iter()
                .any(|read_exon| read_exon.overlap_len(**isoform_intron) > 0);
            inside_read_span || overlaps_read_exon
        })
        .count()
}

fn exon_len_within_span(tx: &Transcript, span: Interval) -> u64 {
    tx.exons
        .iter()
        .map(|exon| u64::from(exon.overlap_len(span)))
        .sum()
}

fn exon_count_within_span(tx: &Transcript, span: Interval) -> usize {
    tx.exons
        .iter()
        .filter(|exon| exon.overlap_len(span) > 0)
        .count()
}

fn read_span(read: &Transcript) -> Interval {
    Interval {
        start: read.tx_start,
        end: read.tx_end,
    }
}

fn assignment_score(
    read: &Transcript,
    isoform: &Transcript,
    junction_offset: u32,
) -> AssignmentScore {
    let read_introns = read.introns();
    let isoform_introns = isoform.introns();
    let (matched_read_introns, matched_isoform_introns) =
        matching_intron_indices(&read_introns, &isoform_introns, junction_offset);

    let read_len = transcript_exon_len(read);
    let isoform_len = exon_len_within_span(isoform, read_span(read));
    let overlap = transcript_exon_overlap(read, isoform);

    AssignmentScore {
        chrom_mismatch: read.chrom != isoform.chrom,
        strand_mismatch: read.strand != isoform.strand,
        unmatched_read_junctions: read_introns
            .len()
            .saturating_sub(matched_read_introns.len()),
        extra_isoform_junctions: covered_unmatched_isoform_junctions(
            read,
            &isoform_introns,
            &matched_isoform_introns,
        ),
        terminal_delta: u64::from(read.tx_start.get().abs_diff(isoform.tx_start.get()))
            + u64::from(read.tx_end.get().abs_diff(isoform.tx_end.get())),
        exon_symmetric_difference: read_len + isoform_len - 2 * overlap,
        exon_count_delta: read
            .exons
            .len()
            .abs_diff(exon_count_within_span(isoform, read_span(read))),
    }
}

fn build_catalog_span_indices(isoforms: &[Transcript]) -> HashMap<CatalogKey, SpanCandidateIndex> {
    let mut grouped: HashMap<CatalogKey, Vec<usize>> = HashMap::new();
    for (idx, isoform) in isoforms.iter().enumerate() {
        grouped
            .entry(CatalogKey::new(None, &isoform.chrom, isoform.strand))
            .or_default()
            .push(idx);
        let mut seen_genes = HashSet::new();
        for gene in transcript_gene_tokens(isoform) {
            if seen_genes.insert(gene) {
                grouped
                    .entry(CatalogKey::new(Some(gene), &isoform.chrom, isoform.strand))
                    .or_default()
                    .push(idx);
            }
        }
    }

    grouped
        .into_iter()
        .map(|(key, indices)| (key, SpanCandidateIndex::build(isoforms, &indices)))
        .collect()
}

fn next_catalog_stamp(seen: &mut [u32], stamp: &mut u32) -> u32 {
    *stamp = stamp.wrapping_add(1);
    if *stamp == 0 {
        seen.fill(0);
        *stamp = 1;
    }
    *stamp
}

fn collect_catalog_candidates_for_key(
    catalog_indices: &HashMap<CatalogKey, SpanCandidateIndex>,
    key: CatalogKey,
    read: &Transcript,
    isoforms: &[Transcript],
    seen: &mut [u32],
    stamp: u32,
    out: &mut Vec<usize>,
) {
    let Some(index) = catalog_indices.get(&key) else {
        return;
    };
    index.collect_overlapping(read, isoforms, seen, stamp, out);
}

pub fn select_unique_best_read_to_isoform(
    reads: &[Transcript],
    isoforms: &[Transcript],
    read_to_isoform: &[(String, String)],
) -> anyhow::Result<Vec<(String, String)>> {
    select_unique_best_read_to_isoform_with_options(
        reads,
        isoforms,
        read_to_isoform,
        UniqueAssignmentOptions::default(),
    )
}

pub fn select_unique_best_read_to_isoform_with_options(
    reads: &[Transcript],
    isoforms: &[Transcript],
    read_to_isoform: &[(String, String)],
    options: UniqueAssignmentOptions,
) -> anyhow::Result<Vec<(String, String)>> {
    let reads_by_name = read_structures_by_name(reads)?;

    let mut isoforms_by_name: HashMap<&str, usize> = HashMap::with_capacity(isoforms.len());
    for (isoform_idx, isoform) in isoforms.iter().enumerate() {
        if let Some(previous_idx) = isoforms_by_name.insert(isoform.name.as_str(), isoform_idx) {
            anyhow::bail!(
                "duplicate isoform id {:?} at catalog indices {previous_idx} and {isoform_idx}; unique assignment requires globally unique isoform ids",
                isoform.name
            );
        }
    }
    let retained_regions = retained_regions_by_isoform(isoforms);
    let catalog_indices = build_catalog_span_indices(isoforms);
    let mut genes_present: HashSet<&str> = HashSet::new();
    for isoform in isoforms {
        for gene in transcript_gene_tokens(isoform) {
            genes_present.insert(gene);
        }
    }

    let mut grouped: BTreeMap<&str, Vec<&str>> = BTreeMap::new();
    let mut seen_pairs: HashSet<(&str, &str)> = HashSet::new();
    for (read_id, isoform_id) in read_to_isoform {
        if !isoforms_by_name.contains_key(isoform_id.as_str()) {
            anyhow::bail!(
                "read_to_isoform references isoform id {isoform_id:?} that is missing from isoform BED"
            );
        }
        if seen_pairs.insert((read_id.as_str(), isoform_id.as_str())) {
            grouped
                .entry(read_id.as_str())
                .or_default()
                .push(isoform_id.as_str());
        }
    }

    // Keep the mapping as the counted-read boundary; catalog expansion should not
    // resurrect unused reads or reads excluded by per-gene downsampling.
    let mapped_read_ids: Vec<&str> = grouped.keys().copied().collect();
    let mut catalog_seen = vec![0u32; isoforms.len()];
    let mut catalog_stamp = 0u32;
    let mut catalog_candidate_indices = Vec::new();
    for read_id in mapped_read_ids {
        let Some(read) = reads_by_name.get(read_id).copied() else {
            continue;
        };
        catalog_candidate_indices.clear();
        let stamp = next_catalog_stamp(&mut catalog_seen, &mut catalog_stamp);
        let read_genes = transcript_gene_tokens(read).collect::<Vec<_>>();
        let read_has_annotated_gene = !read_genes.is_empty();
        let mut domain_found = false;
        if read_has_annotated_gene {
            let mut seen_genes = HashSet::new();
            for gene in read_genes {
                if seen_genes.insert(gene) && genes_present.contains(gene) {
                    domain_found = true;
                    collect_catalog_candidates_for_key(
                        &catalog_indices,
                        CatalogKey::new(Some(gene), &read.chrom, read.strand),
                        read,
                        isoforms,
                        &mut catalog_seen,
                        stamp,
                        &mut catalog_candidate_indices,
                    );
                }
            }
        } else if let Some(mapped_candidates) = grouped.get(read_id) {
            let mut seen_genes: HashSet<&str> = HashSet::new();
            for isoform_id in mapped_candidates {
                let Some(isoform) = isoforms_by_name
                    .get(*isoform_id)
                    .map(|&isoform_idx| &isoforms[isoform_idx])
                else {
                    continue;
                };
                for gene in transcript_gene_tokens(isoform) {
                    if seen_genes.insert(gene) {
                        domain_found = true;
                        collect_catalog_candidates_for_key(
                            &catalog_indices,
                            CatalogKey::new(Some(gene), &read.chrom, read.strand),
                            read,
                            isoforms,
                            &mut catalog_seen,
                            stamp,
                            &mut catalog_candidate_indices,
                        );
                    }
                }
            }
        }

        if !read_has_annotated_gene && !domain_found {
            collect_catalog_candidates_for_key(
                &catalog_indices,
                CatalogKey::new(None, &read.chrom, read.strand),
                read,
                isoforms,
                &mut catalog_seen,
                stamp,
                &mut catalog_candidate_indices,
            );
        }
        if catalog_candidate_indices.is_empty() {
            continue;
        }
        let candidates = grouped
            .get_mut(read_id)
            .expect("read id collected from grouped keys");
        for &isoform_idx in &catalog_candidate_indices {
            let isoform = &isoforms[isoform_idx];
            if catalog_assignment_candidate(read, isoform, options.junction_offset) {
                candidates.push(isoform.name.as_str());
            }
        }
    }

    let mut unique_pairs = Vec::with_capacity(grouped.len());
    for (read_id, mut candidates) in grouped {
        if let Some(read) = reads_by_name.get(read_id).copied() {
            candidates.retain(|isoform_id| {
                isoforms_by_name
                    .get(*isoform_id)
                    .copied()
                    .is_some_and(|isoform_idx| {
                        catalog_assignment_candidate(
                            read,
                            &isoforms[isoform_idx],
                            options.junction_offset,
                        ) && (retained_regions[isoform_idx].is_empty()
                            || retained_regions_supported_by_read(
                                read,
                                &retained_regions[isoform_idx],
                            ))
                    })
            });
        }

        candidates.sort_unstable();
        candidates.dedup();
        if candidates.is_empty() {
            continue;
        }
        if candidates.len() == 1 {
            unique_pairs.push((read_id.to_owned(), candidates[0].to_owned()));
            continue;
        }

        let read = reads_by_name.get(read_id).copied().ok_or_else(|| {
            anyhow::anyhow!(
                "unique assignment needs read structure for multi-mapped read {read_id:?}; \
provide the matching --reads/manifest reads BED"
            )
        })?;

        let mut best: Option<(AssignmentScore, &str)> = None;
        for isoform_id in candidates {
            let isoform_idx = isoforms_by_name
                .get(isoform_id)
                .copied()
                .expect("validated above");
            let score = assignment_score(read, &isoforms[isoform_idx], options.junction_offset);
            if best.as_ref().is_none_or(|(best_score, best_isoform_id)| {
                &score < best_score
                    || (&score == best_score
                        && (!seen_pairs.contains(&(read_id, isoform_id)), isoform_id)
                            < (
                                !seen_pairs.contains(&(read_id, *best_isoform_id)),
                                *best_isoform_id,
                            ))
            }) {
                best = Some((score, isoform_id));
            }
        }
        let (_, best_isoform_id) = best.expect("multi-candidate read has candidates");
        unique_pairs.push((read_id.to_owned(), best_isoform_id.to_owned()));
    }

    Ok(unique_pairs)
}

pub fn count_by_read_to_isoform(
    isoforms: &[Transcript],
    read_to_isoform: &[(String, String)],
) -> anyhow::Result<Vec<CountRecord>> {
    crate::identity::validate_isoform_ids(isoforms)?;
    let isoform_ids: HashSet<&str> = isoforms.iter().map(|tx| tx.name.as_str()).collect();
    let mut unique_pairs: std::collections::BTreeSet<(&str, &str)> =
        std::collections::BTreeSet::new();
    for (read_id, isoform_id) in read_to_isoform {
        if read_id.is_empty() {
            anyhow::bail!("read_to_isoform contains an empty read id");
        }
        if !isoform_ids.contains(isoform_id.as_str()) {
            anyhow::bail!(
                "read_to_isoform references isoform id {isoform_id:?} that is missing from isoform BED"
            );
        }
        unique_pairs.insert((read_id.as_str(), isoform_id.as_str()));
    }

    let mut read_occurrence: HashMap<&str, u32> = HashMap::new();
    for &(read_id, _) in &unique_pairs {
        *read_occurrence.entry(read_id).or_insert(0) += 1;
    }

    let mut counts: HashMap<&str, f64> = HashMap::new();
    for &(read_id, isoform_id) in &unique_pairs {
        let denom = read_occurrence.get(read_id).copied().unwrap_or(0);
        if denom == 0 {
            continue;
        }
        *counts.entry(isoform_id).or_insert(0.0) += 1.0f64 / denom as f64;
    }

    Ok(isoforms
        .iter()
        .map(|isoform| CountRecord {
            gene: gene_name(isoform).unwrap_or("none").to_owned(),
            isoform_id: isoform.name.clone(),
            count: counts.get(isoform.name.as_str()).copied().unwrap_or(0.0),
        })
        .collect())
}

pub fn write_counts_csv<P: AsRef<Path>>(
    path: P,
    records: &[CountRecord],
) -> Result<(), csv::Error> {
    let file = std::fs::File::create(path)?;
    let mut writer = std::io::BufWriter::new(file);
    write_counts_csv_to_writer(&mut writer, records)
}

/// Serialize count records as CSV to an existing writer.
pub fn write_counts_csv_to_writer<W: std::io::Write>(
    writer: &mut W,
    records: &[CountRecord],
) -> Result<(), csv::Error> {
    let mut writer = csv::WriterBuilder::new().from_writer(writer);
    writer.write_record(["gene", "isoform_id", "count"])?;
    for record in records {
        writer.write_record([
            record.gene.as_str(),
            record.isoform_id.as_str(),
            &record.count.to_string(),
        ])?;
    }
    writer.flush()?;
    Ok(())
}

/// Serialize the effective unique-assignment policy as a small, stable TSV.
pub fn unique_assignment_provenance_tsv(options: UniqueAssignmentOptions) -> String {
    format!(
        "format_version\t1\nassignment_mode\tunique\nunique_assignment_junction_offset\t{}\nintron_matcher\tordered_one_to_one_max_cardinality_min_delta\nmicrofeature_collapse\tfalse\n",
        options.junction_offset
    )
}

/// Write the effective unique-assignment policy next to a derived output.
pub fn write_unique_assignment_provenance<P: AsRef<Path>>(
    path: P,
    options: UniqueAssignmentOptions,
) -> Result<(), std::io::Error> {
    std::fs::write(path, unique_assignment_provenance_tsv(options))
}

/// Write effective unique-assignment policy to an existing writer.
pub fn write_unique_assignment_provenance_to_writer<W: std::io::Write>(
    writer: &mut W,
    options: UniqueAssignmentOptions,
) -> Result<(), std::io::Error> {
    writer.write_all(unique_assignment_provenance_tsv(options).as_bytes())
}

#[cfg(test)]
mod tests {
    use crate::model::{Bed12Attrs, Coord, Interval, Strand, Transcript};

    use super::*;

    fn make_tx(name: &str, exons: &[(u32, u32)], name2: &str) -> Transcript {
        let tx_start = exons.iter().map(|(s, _)| *s).min().unwrap_or(0);
        let tx_end = exons.iter().map(|(_, e)| *e).max().unwrap_or(0);
        let exons = exons
            .iter()
            .map(|(s, e)| Interval::new(Coord::new(*s), Coord::new(*e)).unwrap())
            .collect::<Vec<_>>();

        Transcript::new(
            "chr1".to_owned(),
            Strand::Plus,
            Coord::new(tx_start),
            Coord::new(tx_end),
            name.to_owned(),
            exons,
            Bed12Attrs {
                score: 0,
                thick_start: Coord::new(tx_start),
                thick_end: Coord::new(tx_end),
                item_rgb: "0".to_owned(),
                extra_fields: vec![name2.to_owned()],
            },
        )
        .unwrap()
    }

    fn make_tx_with_gene(name: &str, exons: &[(u32, u32)], name2: &str, gene: &str) -> Transcript {
        make_tx_strand_with_gene(name, Strand::Plus, exons, name2, gene)
    }

    fn make_tx_strand_with_gene(
        name: &str,
        strand: Strand,
        exons: &[(u32, u32)],
        name2: &str,
        gene: &str,
    ) -> Transcript {
        make_tx_in_domain(name, "chr1", strand, exons, name2, Some(gene))
    }

    fn make_tx_in_domain(
        name: &str,
        chrom: &str,
        strand: Strand,
        exons: &[(u32, u32)],
        name2: &str,
        gene: Option<&str>,
    ) -> Transcript {
        let tx_start = exons.iter().map(|(s, _)| *s).min().unwrap_or(0);
        let tx_end = exons.iter().map(|(_, e)| *e).max().unwrap_or(0);
        let exons = exons
            .iter()
            .map(|(s, e)| Interval::new(Coord::new(*s), Coord::new(*e)).unwrap())
            .collect::<Vec<_>>();

        let extra_fields = match gene {
            Some(gene) => vec![
                name2.to_owned(),
                "none".to_owned(),
                "none".to_owned(),
                "none".to_owned(),
                "none".to_owned(),
                gene.to_owned(),
            ],
            None => vec![name2.to_owned()],
        };

        Transcript::new(
            chrom.to_owned(),
            strand,
            Coord::new(tx_start),
            Coord::new(tx_end),
            name.to_owned(),
            exons,
            Bed12Attrs {
                score: 0,
                thick_start: Coord::new(tx_start),
                thick_end: Coord::new(tx_end),
                item_rgb: "0".to_owned(),
                extra_fields,
            },
        )
        .unwrap()
    }

    #[test]
    fn counts_split_duplicates_across_isoforms() {
        let references = vec![make_tx("ref1", &[(0, 10)], "ref1")];
        let isoforms = vec![
            make_tx("iso1", &[(0, 10)], "r1,r2,|0"),
            make_tx("iso2", &[(0, 10)], "r2,|0"),
            make_tx("iso3", &[(0, 10)], ",|0"), // empty
        ];

        let records = count_by_subreads(&isoforms, &references).unwrap();
        let iso1 = records.iter().find(|r| r.isoform_id == "iso1").unwrap();
        let iso2 = records.iter().find(|r| r.isoform_id == "iso2").unwrap();
        let iso3 = records.iter().find(|r| r.isoform_id == "iso3").unwrap();

        assert!((iso1.count - 1.5).abs() < 1e-9);
        assert!((iso2.count - 0.5).abs() < 1e-9);
        assert!((iso3.count - 0.0).abs() < 1e-9);
    }

    #[test]
    fn mapping_counts_match_subread_counts() {
        let references = vec![make_tx("ref1", &[(0, 10)], "ref1")];
        let isoforms = vec![
            make_tx("iso1", &[(0, 10)], "r1,r2,|0"),
            make_tx("iso2", &[(0, 10)], "r2,|0"),
        ];
        let pairs = vec![
            ("r1".to_owned(), "iso1".to_owned()),
            ("r2".to_owned(), "iso1".to_owned()),
            ("r2".to_owned(), "iso2".to_owned()),
        ];

        let by_subreads = count_by_subreads(&isoforms, &references).unwrap();
        let by_mapping = count_by_read_to_isoform(&isoforms, &pairs).unwrap();
        assert_eq!(by_subreads.len(), by_mapping.len());

        for (left, right) in by_subreads.iter().zip(by_mapping.iter()) {
            assert_eq!(left.isoform_id, right.isoform_id);
            assert!((left.count - right.count).abs() < 1e-9);
        }
    }

    #[test]
    fn duplicate_mapping_rows_are_idempotent_for_molecule_counts() {
        let isoforms = vec![
            make_tx_with_gene("iso1", &[(0, 10)], "none", "GENEA"),
            make_tx_with_gene("iso2", &[(0, 10)], "none", "GENEA"),
        ];
        let pairs = vec![
            ("molecule".to_owned(), "iso1".to_owned()),
            ("molecule".to_owned(), "iso1".to_owned()),
            ("molecule".to_owned(), "iso2".to_owned()),
        ];

        let counts = count_by_read_to_isoform(&isoforms, &pairs).unwrap();
        assert_eq!(counts[0].count, 0.5);
        assert_eq!(counts[1].count, 0.5);
        assert_eq!(counts.iter().map(|row| row.count).sum::<f64>(), 1.0);
    }

    #[test]
    fn mapping_tsv_round_trips_boundary_whitespace_in_ids() {
        let path = std::env::temp_dir().join(format!(
            "trackcluster-mapping-whitespace-{}-{}.tsv",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let pairs = vec![
            ("r".to_owned(), "iso1".to_owned()),
            (" r ".to_owned(), "iso1".to_owned()),
        ];
        crate::cluster::output::write_read_to_isoform_tsv(&path, &pairs).unwrap();

        let round_tripped = read_read_to_isoform_tsv(&path).unwrap();
        assert_eq!(round_tripped, pairs);
        let isoforms = vec![make_tx_with_gene("iso1", &[(0, 10)], "none", "GENEA")];
        let counts = count_by_read_to_isoform(&isoforms, &round_tripped).unwrap();
        assert_eq!(counts[0].count, 2.0);
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn fractional_count_rejects_duplicate_catalog_ids() {
        let isoforms = vec![
            make_tx("duplicate", &[(0, 10)], "none"),
            make_tx("duplicate", &[(20, 30)], "none"),
        ];
        let error = count_by_read_to_isoform(&isoforms, &[]).unwrap_err();
        assert!(error.to_string().contains("duplicate isoform id"));
    }

    #[test]
    fn count_csv_writes_gene_and_escapes_csv_fields() {
        let path = std::env::temp_dir().join(format!(
            "trackcluster-count-csv-{}-{}.csv",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let rows = vec![CountRecord {
            gene: "GENE,\"quoted\"".to_owned(),
            isoform_id: "iso,one".to_owned(),
            count: 1.25,
        }];

        write_counts_csv(&path, &rows).unwrap();
        let mut reader = csv::Reader::from_path(&path).unwrap();
        assert_eq!(
            reader.headers().unwrap().iter().collect::<Vec<_>>(),
            ["gene", "isoform_id", "count"]
        );
        let record = reader.records().next().unwrap().unwrap();
        assert_eq!(&record[0], "GENE,\"quoted\"");
        assert_eq!(&record[1], "iso,one");
        assert_eq!(&record[2], "1.25");
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn builds_mapping_without_suppressing_read_reference_name_collisions() {
        let references = vec![make_tx("ref1", &[(0, 10)], "ref1")];
        let isoforms = vec![
            make_tx("iso1", &[(0, 10)], "r1,r2,|0"),
            make_tx("iso2", &[(0, 10)], "r2,ref1,|0"),
        ];

        let pairs = read_to_isoform_from_subreads(&isoforms, &references).unwrap();
        assert_eq!(
            pairs,
            vec![
                ("r1".to_owned(), "iso1".to_owned()),
                ("r2".to_owned(), "iso1".to_owned()),
                ("r2".to_owned(), "iso2".to_owned()),
                ("ref1".to_owned(), "iso2".to_owned()),
            ]
        );
    }

    #[test]
    fn unique_assignment_selects_closest_isoform() {
        let reads = vec![make_tx("r1", &[(100, 110), (200, 210)], "none")];
        let isoforms = vec![
            make_tx("long_ref", &[(50, 60), (100, 110), (200, 210)], "none"),
            make_tx("closest_novel", &[(100, 110), (200, 210)], "none"),
        ];
        let pairs = vec![
            ("r1".to_owned(), "long_ref".to_owned()),
            ("r1".to_owned(), "closest_novel".to_owned()),
        ];

        let unique = select_unique_best_read_to_isoform(&reads, &isoforms, &pairs).unwrap();
        assert_eq!(unique, vec![("r1".to_owned(), "closest_novel".to_owned())]);

        let counts = count_by_read_to_isoform(&isoforms, &unique).unwrap();
        let long_ref = counts
            .iter()
            .find(|record| record.isoform_id == "long_ref")
            .unwrap();
        let closest = counts
            .iter()
            .find(|record| record.isoform_id == "closest_novel")
            .unwrap();
        assert_eq!(long_ref.count, 0.0);
        assert_eq!(closest.count, 1.0);
    }

    #[test]
    fn unique_assignment_rejects_conflicting_alignments_for_one_molecule_id() {
        let reads = vec![
            make_tx("r1", &[(100, 110)], "none"),
            make_tx("r1", &[(100, 110), (200, 210)], "none"),
        ];
        let isoforms = vec![
            make_tx("short_isoform", &[(100, 110)], "none"),
            make_tx("long_isoform", &[(100, 110), (200, 210)], "none"),
        ];
        let pairs = vec![
            ("r1".to_owned(), "short_isoform".to_owned()),
            ("r1".to_owned(), "long_isoform".to_owned()),
        ];

        let error = select_unique_best_read_to_isoform(&reads, &isoforms, &pairs).unwrap_err();
        assert!(error.to_string().contains("conflicting alignments"));
        assert!(error.to_string().contains("one unambiguous structure"));
    }

    #[test]
    fn unique_assignment_expands_to_closer_catalog_isoform() {
        let reads = vec![make_tx("r1", &[(100, 110), (200, 210)], "none")];
        let isoforms = vec![
            make_tx("long_ref", &[(50, 60), (100, 110), (200, 210)], "none"),
            make_tx("closest_novel", &[(100, 110), (200, 210)], "none"),
        ];
        let pairs = vec![("r1".to_owned(), "long_ref".to_owned())];

        let unique = select_unique_best_read_to_isoform(&reads, &isoforms, &pairs).unwrap();
        assert_eq!(unique, vec![("r1".to_owned(), "closest_novel".to_owned())]);
    }

    #[test]
    fn unique_assignment_prefers_original_mapping_on_exact_catalog_tie() {
        let reads = vec![make_tx_with_gene(
            "read_b",
            &[(100, 110), (200, 210)],
            "none",
            "GENEA",
        )];
        let isoforms = vec![
            make_tx_with_gene(
                "a_lexical_first",
                &[(100, 110), (200, 210)],
                "none",
                "GENEA",
            ),
            make_tx_with_gene("read_b", &[(100, 110), (200, 210)], "none", "GENEA"),
        ];
        let pairs = vec![("read_b".to_owned(), "read_b".to_owned())];

        let unique = select_unique_best_read_to_isoform(&reads, &isoforms, &pairs).unwrap();
        assert_eq!(unique, vec![("read_b".to_owned(), "read_b".to_owned())]);
    }

    #[test]
    fn unique_assignment_uses_mapped_isoform_gene_when_read_has_no_gene() {
        let reads = vec![make_tx("r1", &[(100, 110), (200, 210)], "none")];
        let isoforms = vec![
            make_tx_with_gene(
                "mapped_long",
                &[(50, 60), (100, 110), (200, 210)],
                "none",
                "GENEA",
            ),
            make_tx_with_gene("z_closest", &[(100, 110), (200, 210)], "none", "GENEA"),
            make_tx_with_gene("a_decoy", &[(100, 110), (200, 210)], "none", "GENEB"),
        ];
        let pairs = vec![("r1".to_owned(), "mapped_long".to_owned())];

        let unique = select_unique_best_read_to_isoform(&reads, &isoforms, &pairs).unwrap();
        assert_eq!(unique, vec![("r1".to_owned(), "z_closest".to_owned())]);
    }

    #[test]
    fn unique_assignment_intersects_multi_gene_metadata_without_cross_gene_leakage() {
        let reads = vec![make_tx_with_gene(
            "multi_gene_read",
            &[(100, 110), (200, 210)],
            "none",
            "GENEA||GENEB",
        )];
        let isoforms = vec![
            make_tx_with_gene(
                "mapped_gene_a",
                &[(50, 60), (100, 110), (200, 210)],
                "none",
                "GENEA",
            ),
            make_tx_with_gene("closest_gene_b", &[(100, 110), (200, 210)], "none", "GENEB"),
            make_tx_with_gene(
                "cross_gene_decoy",
                &[(100, 110), (200, 210)],
                "none",
                "GENEC",
            ),
        ];
        let pairs = vec![("multi_gene_read".to_owned(), "mapped_gene_a".to_owned())];

        assert!(catalog_assignment_candidate(&reads[0], &isoforms[0], 15));
        assert!(catalog_assignment_candidate(&reads[0], &isoforms[1], 15));
        assert!(!catalog_assignment_candidate(&reads[0], &isoforms[2], 15));

        let unique = select_unique_best_read_to_isoform(&reads, &isoforms, &pairs).unwrap();
        assert_eq!(
            unique,
            vec![("multi_gene_read".to_owned(), "closest_gene_b".to_owned())]
        );
    }

    #[test]
    fn unique_assignment_does_not_fallback_to_locus_when_read_gene_exists() {
        let reads = vec![make_tx_with_gene("r1", &[(100, 110)], "none", "GENEA")];
        let isoforms = vec![
            make_tx_with_gene("mapped_far", &[(1000, 1010)], "none", "GENEA"),
            make_tx_with_gene("wrong_gene_overlap", &[(100, 110)], "none", "GENEB"),
        ];
        let pairs = vec![("r1".to_owned(), "mapped_far".to_owned())];

        let unique = select_unique_best_read_to_isoform(&reads, &isoforms, &pairs).unwrap();
        assert!(unique.is_empty());
    }

    #[test]
    fn unique_assignment_prefers_minus_strand_three_prime_early_stop_isoform() {
        let reads = vec![make_tx_strand_with_gene(
            "r_early_stop",
            Strand::Minus,
            &[(160, 200), (300, 350), (400, 500)],
            "none",
            "GENEA",
        )];
        let isoforms = vec![
            make_tx_strand_with_gene(
                "long_ref",
                Strand::Minus,
                &[(100, 200), (300, 350), (400, 500)],
                "none",
                "GENEA",
            ),
            make_tx_strand_with_gene(
                "early_stop",
                Strand::Minus,
                &[(160, 200), (300, 350), (400, 500)],
                "none",
                "GENEA",
            ),
        ];
        let pairs = vec![("r_early_stop".to_owned(), "long_ref".to_owned())];

        let unique = select_unique_best_read_to_isoform(&reads, &isoforms, &pairs).unwrap();
        assert_eq!(
            unique,
            vec![("r_early_stop".to_owned(), "early_stop".to_owned())]
        );
    }

    #[test]
    fn unique_assignment_drops_junctionless_internal_retained_exon_read() {
        let reads = vec![make_tx_with_gene(
            "r_internal",
            &[(140, 160)],
            "none",
            "GENEA",
        )];
        let isoforms = vec![
            make_tx_with_gene(
                "retained",
                &[(100, 110), (120, 200), (300, 310)],
                "none",
                "GENEA",
            ),
            make_tx_with_gene(
                "spliced",
                &[(100, 110), (120, 130), (180, 200), (300, 310)],
                "none",
                "GENEA",
            ),
        ];
        let pairs = vec![("r_internal".to_owned(), "retained".to_owned())];

        let unique = select_unique_best_read_to_isoform(&reads, &isoforms, &pairs).unwrap();
        assert!(unique.is_empty());
    }

    #[test]
    fn unique_assignment_keeps_junctionless_terminal_exon_read() {
        let reads = vec![make_tx_with_gene(
            "r_terminal",
            &[(302, 308)],
            "none",
            "GENEA",
        )];
        let isoforms = vec![make_tx_with_gene(
            "retained",
            &[(100, 110), (120, 200), (300, 310)],
            "none",
            "GENEA",
        )];
        let pairs = vec![("r_terminal".to_owned(), "retained".to_owned())];

        let unique = select_unique_best_read_to_isoform(&reads, &isoforms, &pairs).unwrap();
        assert_eq!(
            unique,
            vec![("r_terminal".to_owned(), "retained".to_owned())]
        );
    }

    #[test]
    fn unique_assignment_keeps_minus_strand_junctionless_three_prime_read() {
        let reads = vec![make_tx_strand_with_gene(
            "r_minus_terminal",
            Strand::Minus,
            &[(102, 108)],
            "none",
            "GENEA",
        )];
        let isoforms = vec![make_tx_strand_with_gene(
            "minus_isoform",
            Strand::Minus,
            &[(100, 110), (200, 210)],
            "none",
            "GENEA",
        )];
        let pairs = vec![("r_minus_terminal".to_owned(), "minus_isoform".to_owned())];

        let unique = select_unique_best_read_to_isoform(&reads, &isoforms, &pairs).unwrap();
        assert_eq!(
            unique,
            vec![("r_minus_terminal".to_owned(), "minus_isoform".to_owned())]
        );
    }

    #[test]
    fn unique_assignment_requires_retained_intron_boundary_support() {
        let reads = vec![make_tx_with_gene(
            "r_suffix",
            &[(300, 310), (400, 410)],
            "none",
            "GENEA",
        )];
        let isoforms = vec![
            make_tx_with_gene(
                "Novel_retained",
                &[(100, 200), (300, 310), (400, 410)],
                "none",
                "GENEA",
            ),
            make_tx_with_gene(
                "Ref_spliced",
                &[(100, 130), (180, 200), (300, 310), (400, 410)],
                "none",
                "GENEA",
            ),
        ];
        let pairs = vec![("r_suffix".to_owned(), "Novel_retained".to_owned())];

        let unique = select_unique_best_read_to_isoform(&reads, &isoforms, &pairs).unwrap();
        assert_eq!(
            unique,
            vec![("r_suffix".to_owned(), "Ref_spliced".to_owned())]
        );
    }

    #[test]
    fn unique_assignment_prefers_retained_isoform_when_read_spans_retained_boundary() {
        let reads = vec![make_tx_with_gene(
            "r_retained",
            &[(120, 200), (300, 310)],
            "none",
            "GENEA",
        )];
        let isoforms = vec![
            make_tx_with_gene("Novel_retained", &[(100, 200), (300, 310)], "none", "GENEA"),
            make_tx_with_gene(
                "Ref_spliced",
                &[(100, 130), (180, 200), (300, 310)],
                "none",
                "GENEA",
            ),
        ];
        let pairs = vec![
            ("r_retained".to_owned(), "Novel_retained".to_owned()),
            ("r_retained".to_owned(), "Ref_spliced".to_owned()),
        ];

        let unique = select_unique_best_read_to_isoform(&reads, &isoforms, &pairs).unwrap();
        assert_eq!(
            unique,
            vec![("r_retained".to_owned(), "Novel_retained".to_owned())]
        );
    }

    #[test]
    fn retained_introns_are_partitioned_by_gene_chromosome_and_strand() {
        let retained = make_tx_in_domain(
            "retained_a",
            "chr1",
            Strand::Plus,
            &[(100, 200), (300, 310)],
            "none",
            Some("GENEA"),
        );
        let cross_domain_suppliers = vec![
            make_tx_in_domain(
                "chr2_supplier",
                "chr2",
                Strand::Plus,
                &[(100, 130), (180, 200)],
                "none",
                Some("GENEA"),
            ),
            make_tx_in_domain(
                "minus_supplier",
                "chr1",
                Strand::Minus,
                &[(100, 130), (180, 200)],
                "none",
                Some("GENEA"),
            ),
            make_tx_in_domain(
                "gene_b_supplier",
                "chr1",
                Strand::Plus,
                &[(100, 130), (180, 200)],
                "none",
                Some("GENEB"),
            ),
        ];

        let mut cross_domain_catalog = vec![retained.clone()];
        cross_domain_catalog.extend(cross_domain_suppliers);
        let regions = retained_regions_by_isoform(&cross_domain_catalog);
        assert!(regions[0].is_empty());

        cross_domain_catalog.push(make_tx_in_domain(
            "same_domain_supplier",
            "chr1",
            Strand::Plus,
            &[(100, 130), (180, 200)],
            "none",
            Some("GENEA"),
        ));
        let regions = retained_regions_by_isoform(&cross_domain_catalog);
        assert_eq!(
            regions[0],
            vec![Interval::new(Coord::new(130), Coord::new(180)).unwrap()]
        );
    }

    #[test]
    fn unannotated_retained_introns_do_not_cross_chromosomes_or_strands() {
        let retained = make_tx_in_domain(
            "unannotated_retained",
            "chr1",
            Strand::Plus,
            &[(100, 200)],
            "none",
            None,
        );
        let mut catalog = vec![
            retained,
            make_tx_in_domain(
                "unannotated_chr2",
                "chr2",
                Strand::Plus,
                &[(100, 130), (180, 200)],
                "none",
                None,
            ),
            make_tx_in_domain(
                "unannotated_minus",
                "chr1",
                Strand::Minus,
                &[(100, 130), (180, 200)],
                "none",
                None,
            ),
        ];
        assert!(retained_regions_by_isoform(&catalog)[0].is_empty());

        catalog.push(make_tx_in_domain(
            "unannotated_same_domain",
            "chr1",
            Strand::Plus,
            &[(100, 130), (180, 200)],
            "none",
            None,
        ));
        assert_eq!(
            retained_regions_by_isoform(&catalog)[0],
            vec![Interval::new(Coord::new(130), Coord::new(180)).unwrap()]
        );
    }

    #[test]
    fn unrelated_chromosome_isoform_cannot_change_unique_assignment() {
        let read = make_tx_in_domain(
            "r_suffix",
            "chr1",
            Strand::Plus,
            &[(300, 310), (400, 410)],
            "none",
            Some("GENEA"),
        );
        let isoform_a = make_tx_in_domain(
            "isoform_a",
            "chr1",
            Strand::Plus,
            &[(100, 200), (300, 310), (400, 410)],
            "none",
            Some("GENEA"),
        );
        let pairs = vec![("r_suffix".to_owned(), "isoform_a".to_owned())];
        let expected = vec![("r_suffix".to_owned(), "isoform_a".to_owned())];

        let base = select_unique_best_read_to_isoform(
            std::slice::from_ref(&read),
            std::slice::from_ref(&isoform_a),
            &pairs,
        )
        .unwrap();
        assert_eq!(base, expected);

        let unrelated = make_tx_in_domain(
            "isoform_b",
            "chr2",
            Strand::Plus,
            &[(100, 130), (180, 200)],
            "none",
            Some("GENEB"),
        );
        let expanded =
            select_unique_best_read_to_isoform(&[read], &[isoform_a, unrelated], &pairs).unwrap();
        assert_eq!(expanded, expected);
    }

    #[test]
    fn unique_assignment_rejects_duplicate_isoform_ids() {
        let isoforms = vec![
            make_tx_with_gene("duplicate", &[(100, 110)], "none", "GENEA"),
            make_tx_with_gene("duplicate", &[(200, 210)], "none", "GENEB"),
        ];

        let error = select_unique_best_read_to_isoform(&[], &isoforms, &[]).unwrap_err();
        assert!(error
            .to_string()
            .contains("duplicate isoform id \"duplicate\""));
        assert!(error.to_string().contains("indices 0 and 1"));
    }

    #[test]
    fn retained_assignment_requires_support_for_every_reached_region() {
        let retained_regions = vec![
            Interval::new(Coord::new(130), Coord::new(180)).unwrap(),
            Interval::new(Coord::new(230), Coord::new(280)).unwrap(),
        ];
        let supports_first_but_splices_second = make_tx(
            "splices_second",
            &[(100, 220), (240, 260), (290, 310)],
            "none",
        );
        assert!(!retained_regions_supported_by_read(
            &supports_first_but_splices_second,
            &retained_regions
        ));

        let reaches_only_supported_first = make_tx("first_only", &[(100, 220)], "none");
        assert!(retained_regions_supported_by_read(
            &reaches_only_supported_first,
            &retained_regions
        ));

        let reaches_neither = make_tx("suffix", &[(300, 310)], "none");
        assert!(!retained_regions_supported_by_read(
            &reaches_neither,
            &retained_regions
        ));
    }

    #[test]
    fn one_isoform_intron_cannot_satisfy_two_read_introns() {
        let read = make_tx(
            "two_read_introns",
            &[(90, 100), (110, 115), (125, 135)],
            "none",
        );
        let isoform = make_tx("one_isoform_intron", &[(90, 107), (117, 135)], "none");

        assert_eq!(
            ordered_one_to_one_intron_matches(&read.introns(), &isoform.introns(), 15).len(),
            1
        );
        assert!(!catalog_assignment_candidate(&read, &isoform, 15));
    }

    #[test]
    fn one_read_intron_cannot_satisfy_two_isoform_introns() {
        let read = make_tx("one_read_intron", &[(90, 107), (117, 135)], "none");
        let isoform = make_tx(
            "two_isoform_introns",
            &[(90, 100), (110, 115), (125, 135)],
            "none",
        );

        assert_eq!(
            ordered_one_to_one_intron_matches(&read.introns(), &isoform.introns(), 15).len(),
            1
        );
        assert!(!catalog_assignment_candidate(&read, &isoform, 15));
    }

    #[test]
    fn ordered_intron_matching_is_strand_independent() {
        for strand in [Strand::Plus, Strand::Minus] {
            let read = make_tx_in_domain(
                "read",
                "chr1",
                strand,
                &[(90, 100), (110, 120), (130, 140)],
                "none",
                Some("GENEA"),
            );
            let isoform = make_tx_in_domain(
                "isoform",
                "chr1",
                strand,
                &[(90, 102), (112, 122), (132, 140)],
                "none",
                Some("GENEA"),
            );

            assert_eq!(
                ordered_one_to_one_intron_matches(&read.introns(), &isoform.introns(), 2),
                vec![(0, 0), (1, 1)]
            );
            assert!(catalog_assignment_candidate(&read, &isoform, 2));
        }
    }

    #[test]
    fn ordered_intron_matching_resolves_near_ties_by_total_delta() {
        let read_introns = vec![
            Interval::new(Coord::new(100), Coord::new(110)).unwrap(),
            Interval::new(Coord::new(112), Coord::new(122)).unwrap(),
        ];
        let exact_tie = vec![Interval::new(Coord::new(106), Coord::new(116)).unwrap()];
        assert_eq!(
            ordered_one_to_one_intron_matches(&read_introns, &exact_tie, 15),
            vec![(0, 0)]
        );

        let closer_to_second = vec![Interval::new(Coord::new(108), Coord::new(118)).unwrap()];
        assert_eq!(
            ordered_one_to_one_intron_matches(&read_introns, &closer_to_second, 15),
            vec![(1, 0)]
        );
    }

    #[test]
    fn ordered_intron_matching_does_not_collapse_repeated_microfeatures() {
        let read_introns = vec![
            Interval::new(Coord::new(100), Coord::new(105)).unwrap(),
            Interval::new(Coord::new(107), Coord::new(112)).unwrap(),
        ];
        let isoform_introns = vec![
            Interval::new(Coord::new(102), Coord::new(107)).unwrap(),
            Interval::new(Coord::new(109), Coord::new(114)).unwrap(),
        ];

        assert_eq!(
            ordered_one_to_one_intron_matches(&read_introns, &isoform_introns, 10),
            vec![(0, 0), (1, 1)]
        );
    }

    #[test]
    fn unique_assignment_junction_tolerance_is_configurable() {
        let reads = vec![make_tx_with_gene(
            "r1",
            &[(90, 100), (110, 130)],
            "none",
            "GENEA",
        )];
        let isoforms = vec![make_tx_with_gene(
            "iso1",
            &[(90, 108), (118, 130)],
            "none",
            "GENEA",
        )];
        let pairs = vec![("r1".to_owned(), "iso1".to_owned())];

        let strict = select_unique_best_read_to_isoform_with_options(
            &reads,
            &isoforms,
            &pairs,
            UniqueAssignmentOptions { junction_offset: 7 },
        )
        .unwrap();
        assert!(strict.is_empty());

        let tolerant = select_unique_best_read_to_isoform_with_options(
            &reads,
            &isoforms,
            &pairs,
            UniqueAssignmentOptions { junction_offset: 8 },
        )
        .unwrap();
        assert_eq!(tolerant, pairs);
        assert_eq!(
            UniqueAssignmentOptions::default().junction_offset,
            DEFAULT_UNIQUE_ASSIGNMENT_JUNCTION_OFFSET
        );
        assert!(
            unique_assignment_provenance_tsv(UniqueAssignmentOptions { junction_offset: 8 })
                .contains("unique_assignment_junction_offset\t8\n")
        );
    }

    #[test]
    fn unique_assignment_does_not_create_catalog_only_reads() {
        let reads = vec![make_tx("r1", &[(100, 110), (200, 210)], "none")];
        let isoforms = vec![make_tx("closest_novel", &[(100, 110), (200, 210)], "none")];
        let pairs = Vec::new();

        let unique = select_unique_best_read_to_isoform(&reads, &isoforms, &pairs).unwrap();
        assert!(unique.is_empty());
    }
}
