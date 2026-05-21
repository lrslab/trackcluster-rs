use std::collections::{BTreeMap, HashMap, HashSet};
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

use crate::model::{Interval, Strand, Transcript};

pub mod multi;

const GENE_NAME_COL: usize = 5;
const UNIQUE_ASSIGNMENT_JUNCTION_OFFSET: u32 = 15;

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

fn matched_junction_counts(
    read_introns: &[Interval],
    isoform_introns: &[Interval],
    offset: u32,
) -> (usize, usize) {
    let matched_read_junctions = read_introns
        .iter()
        .filter(|read_intron| {
            isoform_introns
                .iter()
                .any(|isoform_intron| junctions_match(**read_intron, *isoform_intron, offset))
        })
        .count();
    let matched_isoform_junctions = isoform_introns
        .iter()
        .filter(|isoform_intron| {
            read_introns
                .iter()
                .any(|read_intron| junctions_match(*read_intron, **isoform_intron, offset))
        })
        .count();
    (matched_read_junctions, matched_isoform_junctions)
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
        return true;
    }
    if isoform_introns.is_empty() {
        return false;
    }

    read_introns.iter().all(|read_intron| {
        isoform_introns.iter().any(|isoform_intron| {
            junctions_match(
                *read_intron,
                *isoform_intron,
                UNIQUE_ASSIGNMENT_JUNCTION_OFFSET,
            )
        })
    })
}

fn assignment_score(read: &Transcript, isoform: &Transcript) -> AssignmentScore {
    let read_introns = read.introns();
    let isoform_introns = isoform.introns();
    let (matched_read_junctions, matched_isoform_junctions) = matched_junction_counts(
        &read_introns,
        &isoform_introns,
        UNIQUE_ASSIGNMENT_JUNCTION_OFFSET,
    );

    let read_len = transcript_exon_len(read);
    let isoform_len = transcript_exon_len(isoform);
    let overlap = transcript_exon_overlap(read, isoform);

    AssignmentScore {
        chrom_mismatch: read.chrom != isoform.chrom,
        strand_mismatch: read.strand != isoform.strand,
        unmatched_read_junctions: read_introns.len().saturating_sub(matched_read_junctions),
        extra_isoform_junctions: isoform_introns
            .len()
            .saturating_sub(matched_isoform_junctions),
        terminal_delta: u64::from(read.tx_start.get().abs_diff(isoform.tx_start.get()))
            + u64::from(read.tx_end.get().abs_diff(isoform.tx_end.get())),
        exon_symmetric_difference: read_len + isoform_len - 2 * overlap,
        exon_count_delta: read.exons.len().abs_diff(isoform.exons.len()),
    }
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
    let mut isoforms_by_gene: HashMap<&str, Vec<&str>> = HashMap::new();
    let mut isoforms_by_locus: HashMap<(&str, Strand), Vec<&str>> = HashMap::new();
    for isoform in isoforms {
        if let Some(gene) = gene_name(isoform) {
            isoforms_by_gene
                .entry(gene)
                .or_default()
                .push(isoform.name.as_str());
        }
        isoforms_by_locus
            .entry((isoform.chrom.as_str(), isoform.strand))
            .or_default()
            .push(isoform.name.as_str());
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

    for read in reads_by_name.values().copied() {
        let catalog_candidates = gene_name(read)
            .and_then(|gene| isoforms_by_gene.get(gene))
            .or_else(|| isoforms_by_locus.get(&(read.chrom.as_str(), read.strand)));
        let Some(catalog_candidates) = catalog_candidates else {
            continue;
        };
        let candidates = grouped.entry(read.name.as_str()).or_default();
        for isoform_id in catalog_candidates {
            let isoform = isoforms_by_name
                .get(*isoform_id)
                .copied()
                .expect("catalog index derives from isoforms");
            if catalog_assignment_candidate(read, isoform) {
                candidates.push(*isoform_id);
            }
        }
    }

    let mut unique_pairs = Vec::with_capacity(grouped.len());
    for (read_id, mut candidates) in grouped {
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
                (&score, isoform_id) < (best_score, *best_isoform_id)
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
    fn unique_assignment_can_assign_catalog_only_read() {
        let reads = vec![make_tx("r1", &[(100, 110), (200, 210)], "none")];
        let isoforms = vec![make_tx("closest_novel", &[(100, 110), (200, 210)], "none")];
        let pairs = Vec::new();

        let unique = select_unique_best_read_to_isoform(&reads, &isoforms, &pairs).unwrap();
        assert_eq!(unique, vec![("r1".to_owned(), "closest_novel".to_owned())]);
    }
}
