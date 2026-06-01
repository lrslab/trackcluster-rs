use std::collections::{BTreeMap, HashMap, HashSet};
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

use crate::model::{Coord, Interval, Strand, Transcript};

pub mod multi;

const GENE_NAME_COL: usize = 5;
const UNIQUE_ASSIGNMENT_JUNCTION_OFFSET: u32 = 15;
const UNIQUE_CATALOG_BIN_SIZE: u32 = 16_384;

#[derive(Debug, Default, Clone, PartialEq)]
pub struct CountRecord {
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

pub(crate) fn parse_subreads(tx: &Transcript) -> Vec<&str> {
    let Some(name2) = tx.extra_fields.first() else {
        return Vec::new();
    };
    if !name2.contains(',') {
        return Vec::new();
    }
    let sub_part = name2.split(",|").next().unwrap_or(name2.as_str());
    sub_part
        .split(',')
        .map(str::trim)
        .filter(|token| !token.is_empty())
        .collect()
}

pub fn count_by_subreads(isoforms: &[Transcript], references: &[Transcript]) -> Vec<CountRecord> {
    let ref_names: HashSet<&str> = references.iter().map(|tx| tx.name.as_str()).collect();

    let mut occ: HashMap<&str, u32> = HashMap::new();
    for isoform in isoforms {
        for name in parse_subreads(isoform) {
            if ref_names.contains(name) {
                continue;
            }
            *occ.entry(name).or_insert(0) += 1;
        }
    }

    isoforms
        .iter()
        .map(|isoform| {
            let mut coverage = 0.0f64;
            let mut subreads = parse_subreads(isoform);
            subreads.sort_unstable();
            for name in subreads {
                if ref_names.contains(name) {
                    continue;
                }
                let denom = occ.get(name).copied().unwrap_or(0);
                if denom > 0 {
                    coverage += 1.0f64 / denom as f64;
                }
            }

            CountRecord {
                isoform_id: isoform.name.clone(),
                count: coverage,
            }
        })
        .collect()
}

pub fn read_to_isoform_from_subreads(
    isoforms: &[Transcript],
    references: &[Transcript],
) -> Vec<(String, String)> {
    let ref_names: HashSet<&str> = references.iter().map(|tx| tx.name.as_str()).collect();
    let mut pairs = Vec::new();
    for isoform in isoforms {
        for read_name in parse_subreads(isoform) {
            if ref_names.contains(read_name) {
                continue;
            }
            pairs.push((read_name.to_owned(), isoform.name.clone()));
        }
    }
    pairs
}

pub fn read_read_to_isoform_tsv<P: AsRef<Path>>(
    path: P,
) -> Result<Vec<(String, String)>, std::io::Error> {
    let file = std::fs::File::open(path)?;
    let reader = BufReader::new(file);

    let mut pairs: Vec<(String, String)> = Vec::new();
    for (line_no, line) in reader.lines().enumerate() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        let Some((read_id, isoform_id)) = line.split_once('\t') else {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("invalid read_to_isoform line {}: {:?}", line_no + 1, line),
            ));
        };

        let read_id = read_id.trim();
        let isoform_id = isoform_id.trim();
        if read_id.is_empty() || isoform_id.is_empty() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "empty read/isoform value at read_to_isoform line {}",
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

fn read_structures_by_name(reads: &[Transcript]) -> HashMap<&str, &Transcript> {
    let mut reads_by_name: HashMap<&str, &Transcript> = HashMap::new();
    for read in reads {
        let read_id = read.name.as_str();
        match reads_by_name.get(read_id).copied() {
            Some(existing) if transcript_exon_len(read) <= transcript_exon_len(existing) => {}
            _ => {
                reads_by_name.insert(read_id, read);
            }
        }
    }
    reads_by_name
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
    let gene = tx.extra_fields.get(GENE_NAME_COL)?.trim();
    if gene.is_empty() || gene == "none" {
        None
    } else {
        Some(gene)
    }
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

fn retained_regions_by_isoform(isoforms: &[Transcript]) -> HashMap<String, Vec<Interval>> {
    let mut out = HashMap::new();
    let mut gene_introns = Vec::new();
    for isoform in isoforms {
        gene_introns.extend(isoform.introns());
    }
    gene_introns.sort_unstable();
    gene_introns.dedup();

    if gene_introns.is_empty() {
        return out;
    }

    for isoform in isoforms {
        let mut retained_regions = Vec::new();
        for exon in &isoform.exons {
            let mut idx = lower_bound_introns_by_start(&gene_introns, exon.start);
            while idx < gene_introns.len() && gene_introns[idx].start <= exon.end {
                let intron = gene_introns[idx];
                if !intron.is_empty() && interval_contains(*exon, intron) {
                    retained_regions.push(intron);
                }
                idx += 1;
            }
        }

        retained_regions.sort_unstable();
        retained_regions.dedup();
        if !retained_regions.is_empty() {
            out.insert(isoform.name.clone(), retained_regions);
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
    retained_regions
        .iter()
        .any(|region| read_supports_retained_region(read, *region))
}

fn matching_intron_indices(
    read_introns: &[Interval],
    isoform_introns: &[Interval],
    offset: u32,
) -> (HashSet<usize>, HashSet<usize>) {
    let mut matched_read = HashSet::new();
    let mut matched_isoform = HashSet::new();
    for (read_idx, read_intron) in read_introns.iter().enumerate() {
        for (isoform_idx, isoform_intron) in isoform_introns.iter().enumerate() {
            if junctions_match(*read_intron, *isoform_intron, offset) {
                matched_read.insert(read_idx);
                matched_isoform.insert(isoform_idx);
            }
        }
    }
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

fn catalog_assignment_candidate(read: &Transcript, isoform: &Transcript) -> bool {
    if read.chrom != isoform.chrom || read.strand != isoform.strand {
        return false;
    }
    if let (Some(read_gene), Some(isoform_gene)) = (gene_name(read), gene_name(isoform)) {
        if read_gene != isoform_gene {
            return false;
        }
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

    let (matched_read_introns, matched_isoform_introns) = matching_intron_indices(
        &read_introns,
        &isoform_introns,
        UNIQUE_ASSIGNMENT_JUNCTION_OFFSET,
    );

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

fn assignment_score(read: &Transcript, isoform: &Transcript) -> AssignmentScore {
    let read_introns = read.introns();
    let isoform_introns = isoform.introns();
    let (matched_read_introns, matched_isoform_introns) = matching_intron_indices(
        &read_introns,
        &isoform_introns,
        UNIQUE_ASSIGNMENT_JUNCTION_OFFSET,
    );

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
        if let Some(gene) = gene_name(isoform) {
            grouped
                .entry(CatalogKey::new(Some(gene), &isoform.chrom, isoform.strand))
                .or_default()
                .push(idx);
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
    let reads_by_name = read_structures_by_name(reads);

    let isoforms_by_name: HashMap<&str, &Transcript> = isoforms
        .iter()
        .map(|isoform| (isoform.name.as_str(), isoform))
        .collect();
    let retained_regions = retained_regions_by_isoform(isoforms);
    let catalog_indices = build_catalog_span_indices(isoforms);
    let mut genes_present: HashSet<&str> = HashSet::new();
    for isoform in isoforms {
        if let Some(gene) = gene_name(isoform) {
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
        let mut domain_found = false;
        if let Some(gene) = gene_name(read) {
            if genes_present.contains(gene) {
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
        } else if let Some(mapped_candidates) = grouped.get(read_id) {
            let mut seen_genes: HashSet<&str> = HashSet::new();
            for isoform_id in mapped_candidates {
                let Some(gene) = isoforms_by_name
                    .get(*isoform_id)
                    .copied()
                    .and_then(gene_name)
                else {
                    continue;
                };
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

        if !domain_found {
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
            if catalog_assignment_candidate(read, isoform) {
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
                    .is_some_and(|isoform| {
                        catalog_assignment_candidate(read, isoform)
                            && retained_regions.get(*isoform_id).is_none_or(|regions| {
                                retained_regions_supported_by_read(read, regions)
                            })
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
            let isoform = isoforms_by_name
                .get(isoform_id)
                .copied()
                .expect("validated above");
            let score = assignment_score(read, isoform);
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
) -> Vec<CountRecord> {
    let mut read_occurrence: HashMap<&str, u32> = HashMap::new();
    for (read_id, _) in read_to_isoform {
        *read_occurrence.entry(read_id.as_str()).or_insert(0) += 1;
    }

    let mut counts: HashMap<&str, f64> = HashMap::new();
    for (read_id, isoform_id) in read_to_isoform {
        let denom = read_occurrence.get(read_id.as_str()).copied().unwrap_or(0);
        if denom == 0 {
            continue;
        }
        *counts.entry(isoform_id.as_str()).or_insert(0.0) += 1.0f64 / denom as f64;
    }

    isoforms
        .iter()
        .map(|isoform| CountRecord {
            isoform_id: isoform.name.clone(),
            count: counts.get(isoform.name.as_str()).copied().unwrap_or(0.0),
        })
        .collect()
}

pub fn write_counts_csv<P: AsRef<Path>>(
    path: P,
    records: &[CountRecord],
) -> Result<(), std::io::Error> {
    let mut writer = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(&mut writer, "isoform_id,count")?;
    for record in records {
        writeln!(&mut writer, "{},{}", record.isoform_id, record.count)?;
    }
    Ok(())
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
        let tx_start = exons.iter().map(|(s, _)| *s).min().unwrap_or(0);
        let tx_end = exons.iter().map(|(_, e)| *e).max().unwrap_or(0);
        let exons = exons
            .iter()
            .map(|(s, e)| Interval::new(Coord::new(*s), Coord::new(*e)).unwrap())
            .collect::<Vec<_>>();

        Transcript::new(
            "chr1".to_owned(),
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
                extra_fields: vec![
                    name2.to_owned(),
                    "none".to_owned(),
                    "none".to_owned(),
                    "none".to_owned(),
                    "none".to_owned(),
                    gene.to_owned(),
                ],
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

        let records = count_by_subreads(&isoforms, &references);
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

        let by_subreads = count_by_subreads(&isoforms, &references);
        let by_mapping = count_by_read_to_isoform(&isoforms, &pairs);
        assert_eq!(by_subreads.len(), by_mapping.len());

        for (left, right) in by_subreads.iter().zip(by_mapping.iter()) {
            assert_eq!(left.isoform_id, right.isoform_id);
            assert!((left.count - right.count).abs() < 1e-9);
        }
    }

    #[test]
    fn builds_mapping_from_embedded_subreads() {
        let references = vec![make_tx("ref1", &[(0, 10)], "ref1")];
        let isoforms = vec![
            make_tx("iso1", &[(0, 10)], "r1,r2,|0"),
            make_tx("iso2", &[(0, 10)], "r2,ref1,|0"),
        ];

        let pairs = read_to_isoform_from_subreads(&isoforms, &references);
        assert_eq!(
            pairs,
            vec![
                ("r1".to_owned(), "iso1".to_owned()),
                ("r2".to_owned(), "iso1".to_owned()),
                ("r2".to_owned(), "iso2".to_owned()),
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

        let counts = count_by_read_to_isoform(&isoforms, &unique);
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
    fn unique_assignment_uses_longest_duplicate_read_structure() {
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

        let unique = select_unique_best_read_to_isoform(&reads, &isoforms, &pairs).unwrap();
        assert_eq!(unique, vec![("r1".to_owned(), "long_isoform".to_owned())]);
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
    fn unique_assignment_does_not_create_catalog_only_reads() {
        let reads = vec![make_tx("r1", &[(100, 110), (200, 210)], "none")];
        let isoforms = vec![make_tx("closest_novel", &[(100, 110), (200, 210)], "none")];
        let pairs = Vec::new();

        let unique = select_unique_best_read_to_isoform(&reads, &isoforms, &pairs).unwrap();
        assert!(unique.is_empty());
    }
}
