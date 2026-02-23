use std::collections::{HashMap, HashSet};

use crate::cluster::result::ClusterResult;
use crate::interval::{cluster_by_span, exonic_overlap_bp, sort_by_coord, StrandMode};
use crate::model::{Interval, Strand, Transcript};

#[derive(Clone, Debug)]
struct Track {
    tx: Transcript,
    subreads: HashSet<String>,
    exon_len: u32,
    introns: Vec<Interval>,
    intron_len: u32,
}

impl Track {
    fn new(tx: Transcript) -> Self {
        let mut subreads: HashSet<String> = HashSet::new();
        if !is_isoform_anno(&tx) {
            subreads.insert(tx.name.clone());
        }
        let exon_len = tx.exons.iter().map(|exon| exon.len()).sum();
        let introns = tx.introns();
        let intron_len = introns.iter().map(|intron| intron.len()).sum();
        Self {
            tx,
            subreads,
            exon_len,
            introns,
            intron_len,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum DistanceMode {
    Ratio,
    RatioShort,
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct PartitionKey {
    chrom: String,
    strand: Strand,
}

const NAME2_COL: usize = 0;
const TTYPE_COL: usize = 4;

fn set_extra(tx: &mut Transcript, idx: usize, value: String) {
    if tx.extra_fields.len() <= idx {
        tx.extra_fields.resize(idx + 1, "none".to_owned());
    }
    tx.extra_fields[idx] = value;
}

fn ttype(tx: &Transcript) -> Option<&str> {
    tx.extra_fields.get(TTYPE_COL).map(|value| value.as_str())
}

fn is_isoform_anno(tx: &Transcript) -> bool {
    matches!(ttype(tx), Some("isoform_anno"))
}

fn overlap_bp(a: &[Interval], b: &[Interval]) -> u32 {
    let mut total: u32 = 0;
    let mut ai: usize = 0;
    let mut bi: usize = 0;

    while ai < a.len() && bi < b.len() {
        let a_iv = a[ai];
        let b_iv = b[bi];

        if a_iv.end <= b_iv.start {
            ai += 1;
            continue;
        }
        if b_iv.end <= a_iv.start {
            bi += 1;
            continue;
        }

        total += a_iv.overlap_len(b_iv);

        if a_iv.end <= b_iv.end {
            ai += 1;
        } else {
            bi += 1;
        }
    }

    total
}

fn combined_distance(a: &Track, b: &Track, mode: DistanceMode, intron_weight: f64) -> f64 {
    let exon_overlap = exonic_overlap_bp(&a.tx, &b.tx) as f64;
    let a_exon_len = a.exon_len as f64;
    let b_exon_len = b.exon_len as f64;

    let exon_denom = match mode {
        DistanceMode::Ratio => (a_exon_len + b_exon_len - exon_overlap).max(1.0),
        DistanceMode::RatioShort => a_exon_len.min(b_exon_len).max(1.0),
    };
    let exon_dist = 1.0 - (exon_overlap / exon_denom);

    let intron_overlap = overlap_bp(&a.introns, &b.introns) as f64;
    let a_intron_len = a.intron_len as f64;
    let b_intron_len = b.intron_len as f64;

    let intron_dist = if a_intron_len.min(b_intron_len) <= 0.0 {
        0.0
    } else {
        let denom = match mode {
            DistanceMode::Ratio => (a_intron_len + b_intron_len - intron_overlap).max(1.0),
            DistanceMode::RatioShort => a_intron_len.min(b_intron_len).max(1.0),
        };
        1.0 - (intron_overlap / denom)
    };

    (exon_dist + intron_weight * intron_dist) / (1.0 + intron_weight)
}

fn merge_subreads(src: usize, dst: usize, tracks: &mut [Track]) {
    let name = tracks[src].tx.name.clone();
    tracks[dst].subreads.insert(name);

    let subs: Vec<String> = tracks[src].subreads.iter().cloned().collect();
    tracks[dst].subreads.extend(subs);
}

fn readall(tracks: &[Track]) -> HashSet<String> {
    let mut names: HashSet<String> = HashSet::new();
    for track in tracks {
        names.insert(track.tx.name.clone());
        for sub in &track.subreads {
            names.insert(sub.clone());
        }
    }
    names
}

fn readall_subset(tracks: &[Track], keep: &HashSet<usize>) -> HashSet<String> {
    let mut names: HashSet<String> = HashSet::new();
    for &idx in keep {
        let track = &tracks[idx];
        names.insert(track.tx.name.clone());
        for sub in &track.subreads {
            names.insert(sub.clone());
        }
    }
    names
}

fn filter_pass(
    tracks: Vec<Track>,
    mode: DistanceMode,
    cutoff: f64,
    intron_weight: f64,
    sw_score: u32,
) -> Vec<Track> {
    let mut tracks = tracks;
    let mut drop: HashSet<usize> = HashSet::new();

    for i in 0..tracks.len() {
        for j in (i + 1)..tracks.len() {
            let distance = combined_distance(&tracks[i], &tracks[j], mode, intron_weight);
            if distance >= cutoff {
                continue;
            }

            let li = tracks[i].exon_len;
            let lj = tracks[j].exon_len;
            match li.cmp(&lj) {
                std::cmp::Ordering::Less => {
                    if mode == DistanceMode::Ratio || tracks[i].tx.score <= sw_score {
                        drop.insert(i);
                        merge_subreads(i, j, &mut tracks);
                    }
                }
                std::cmp::Ordering::Equal => {
                    drop.insert(i);
                    merge_subreads(i, j, &mut tracks);
                }
                std::cmp::Ordering::Greater => {
                    if mode == DistanceMode::Ratio || tracks[j].tx.score <= sw_score {
                        drop.insert(j);
                        merge_subreads(j, i, &mut tracks);
                    }
                }
            }
        }
    }

    let mut keep: HashSet<usize> = (0..tracks.len()).collect();
    for idx in &drop {
        keep.remove(idx);
    }

    for (idx, track) in tracks.iter().enumerate() {
        if is_isoform_anno(&track.tx) {
            keep.insert(idx);
        }
    }

    let all_names = readall(&tracks);
    let kept_names = readall_subset(&tracks, &keep);
    let missed = all_names
        .difference(&kept_names)
        .cloned()
        .collect::<Vec<_>>();

    if !missed.is_empty() {
        let mut pos: HashMap<String, usize> = HashMap::new();
        for (idx, track) in tracks.iter().enumerate() {
            pos.insert(track.tx.name.clone(), idx);
        }
        for name in missed {
            if let Some(idx) = pos.get(&name).copied() {
                keep.insert(idx);
            }
        }
    }

    let mut keep_vec: Vec<usize> = keep.into_iter().collect();
    keep_vec.sort_unstable();

    keep_vec
        .into_iter()
        .map(|idx| tracks[idx].clone())
        .collect()
}

fn build_read_to_isoform(isoforms: &[Track], ref_names: &HashSet<String>) -> Vec<(String, String)> {
    let mut pairs: Vec<(String, String)> = Vec::new();
    for track in isoforms {
        for subread in &track.subreads {
            if !ref_names.contains(subread) {
                pairs.push((subread.clone(), track.tx.name.clone()));
            }
        }
    }
    pairs.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));
    pairs
}

fn update_name2_with_coverage(isoforms: &mut [Track], ref_names: &HashSet<String>) {
    let values: Vec<String> = {
        let mut occurrence: HashMap<&str, u32> = HashMap::new();
        for track in isoforms.iter() {
            for name in &track.subreads {
                if !ref_names.contains(name) {
                    *occurrence.entry(name.as_str()).or_insert(0) += 1;
                }
            }
        }

        isoforms
            .iter()
            .map(|track| {
                let mut subreads: Vec<&str> = track.subreads.iter().map(|s| s.as_str()).collect();
                subreads.sort_unstable();
                let joined = subreads.join(",");

                let mut coverage = 0.0f64;
                for name in &subreads {
                    if ref_names.contains(*name) {
                        continue;
                    }
                    let denom = occurrence.get(*name).copied().unwrap_or(0);
                    if denom > 0 {
                        coverage += 1.0f64 / denom as f64;
                    }
                }

                format!("{joined},|{coverage}")
            })
            .collect()
    };

    for (track, value) in isoforms.iter_mut().zip(values) {
        set_extra(&mut track.tx, NAME2_COL, value);
    }
}

struct PartitionResult {
    isoforms: Vec<Transcript>,
    pairs: Vec<(String, String)>,
}

struct WorkItem {
    index: usize,
    ref_indices: Vec<usize>,
    read_indices: Vec<usize>,
}

fn process_partition(
    references: &[Transcript],
    reads: &[Transcript],
    ref_indices: &[usize],
    read_indices: &[usize],
    ref_names: &HashSet<String>,
) -> PartitionResult {
    let mut records: Vec<Transcript> = Vec::with_capacity(ref_indices.len() + read_indices.len());
    for &idx in ref_indices {
        records.push(references[idx].clone());
    }
    for &idx in read_indices {
        records.push(reads[idx].clone());
    }

    sort_by_coord(&mut records);
    let loci = cluster_by_span(&records, StrandMode::Match);
    let mut records: Vec<Option<Transcript>> = records.into_iter().map(Some).collect();

    let mut isoforms: Vec<Transcript> = Vec::new();
    let mut pairs: Vec<(String, String)> = Vec::new();

    for locus in loci {
        let mut tracks: Vec<Track> = Vec::with_capacity(locus.members.len());
        for idx in locus.members {
            let tx = records[idx].take().expect("record already consumed");
            tracks.push(Track::new(tx));
        }

        let tracks = filter_pass(tracks, DistanceMode::Ratio, 0.05, 0.5, 11);
        let mut tracks = filter_pass(tracks, DistanceMode::RatioShort, 0.01, 0.5, 11);

        update_name2_with_coverage(&mut tracks, ref_names);
        pairs.extend(build_read_to_isoform(&tracks, ref_names));
        isoforms.extend(tracks.into_iter().map(|track| track.tx));
    }

    PartitionResult { isoforms, pairs }
}

pub fn cluster(
    reads: &[Transcript],
    references: Option<&[Transcript]>,
    threads: usize,
) -> ClusterResult {
    let references = match references {
        Some(references) => references,
        None => {
            return ClusterResult {
                isoforms: Vec::new(),
                read_to_isoform: Vec::new(),
                unused: Vec::new(),
            }
        }
    };

    let threads = threads.max(1);
    let ref_names: std::sync::Arc<HashSet<String>> =
        std::sync::Arc::new(references.iter().map(|tx| tx.name.clone()).collect());

    let mut refs_by_key: HashMap<PartitionKey, Vec<usize>> = HashMap::new();
    for (idx, tx) in references.iter().enumerate() {
        refs_by_key
            .entry(PartitionKey {
                chrom: tx.chrom.clone(),
                strand: tx.strand,
            })
            .or_default()
            .push(idx);
    }

    let mut reads_by_key: HashMap<PartitionKey, Vec<usize>> = HashMap::new();
    for (idx, read) in reads.iter().enumerate() {
        let key = PartitionKey {
            chrom: read.chrom.clone(),
            strand: read.strand,
        };
        if !refs_by_key.contains_key(&key) {
            continue;
        }
        reads_by_key.entry(key).or_default().push(idx);
    }

    let mut all_isoforms: Vec<Transcript> = Vec::new();
    let mut all_pairs: Vec<(String, String)> = Vec::new();

    let mut keys: Vec<PartitionKey> = refs_by_key.keys().cloned().collect();
    keys.sort_by(|a, b| a.chrom.cmp(&b.chrom).then_with(|| a.strand.cmp(&b.strand)));

    let mut work: Vec<WorkItem> = Vec::with_capacity(keys.len());
    for (index, key) in keys.iter().enumerate() {
        work.push(WorkItem {
            index,
            ref_indices: refs_by_key.remove(key).unwrap_or_default(),
            read_indices: reads_by_key.remove(key).unwrap_or_default(),
        });
    }

    let mut parts: Vec<Option<PartitionResult>> = (0..keys.len()).map(|_| None).collect();
    if threads == 1 || work.len() <= 1 {
        for item in work {
            parts[item.index] = Some(process_partition(
                references,
                reads,
                &item.ref_indices,
                &item.read_indices,
                &ref_names,
            ));
        }
    } else {
        use std::sync::{mpsc, Arc, Mutex};

        let queue = Arc::new(Mutex::new(work));
        let (tx, rx) = mpsc::channel::<(usize, PartitionResult)>();

        let worker_count = threads.min(keys.len());
        std::thread::scope(|scope| {
            for _ in 0..worker_count {
                let queue = Arc::clone(&queue);
                let tx = tx.clone();
                let ref_names = Arc::clone(&ref_names);

                scope.spawn(move || loop {
                    let item = {
                        let mut guard = queue.lock().expect("work queue poisoned");
                        guard.pop()
                    };
                    let Some(item) = item else {
                        break;
                    };

                    let result = process_partition(
                        references,
                        reads,
                        &item.ref_indices,
                        &item.read_indices,
                        &ref_names,
                    );
                    if tx.send((item.index, result)).is_err() {
                        break;
                    }
                });
            }
            drop(tx);

            for _ in 0..keys.len() {
                let (index, result) = rx.recv().expect("worker hung up");
                parts[index] = Some(result);
            }
        });
    }

    for part in parts.into_iter().flatten() {
        all_isoforms.extend(part.isoforms);
        all_pairs.extend(part.pairs);
    }

    all_pairs.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));

    ClusterResult {
        isoforms: all_isoforms,
        read_to_isoform: all_pairs,
        unused: Vec::new(),
    }
}

#[cfg(test)]
mod tests {
    use crate::model::{Bed12Attrs, Coord, Interval, Strand, Transcript};

    use super::*;

    fn make_tx(
        name: &str,
        strand: Strand,
        exons: &[(u32, u32)],
        ttype: &str,
        score: u32,
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
                score,
                thick_start: Coord::new(0),
                thick_end: Coord::new(0),
                item_rgb: "0".to_owned(),
                extra_fields: vec![
                    "none".to_owned(),
                    "none".to_owned(),
                    "none".to_owned(),
                    "-1,".to_owned(),
                    ttype.to_owned(),
                    "none".to_owned(),
                    "none".to_owned(),
                    "none".to_owned(),
                ],
            },
        )
        .unwrap()
    }

    #[test]
    fn merges_near_identical_reads_into_longer() {
        let refs = vec![make_tx(
            "ref",
            Strand::Plus,
            &[(100, 110), (120, 130), (140, 150)],
            "isoform_anno",
            100,
        )];
        let reads = vec![make_tx(
            "read",
            Strand::Plus,
            &[(100, 110), (120, 130), (140, 150)],
            "nanopore_read",
            0,
        )];

        let result = cluster(&reads, Some(&refs), 1);
        assert!(!result.isoforms.is_empty());
    }
}
