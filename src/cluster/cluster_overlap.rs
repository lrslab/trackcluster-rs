use std::collections::{HashMap, HashSet};

use crate::cluster::{clusterj::Name2Mode, result::ClusterResult};
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

pub const DEFAULT_CUTOFF1: f64 = 0.05;
pub const DEFAULT_CUTOFF2: f64 = 0.01;
pub const DEFAULT_INTRON_WEIGHT: f64 = 0.5;
pub const DEFAULT_SW_SCORE: i64 = 11;

#[derive(Clone, Copy, Debug)]
pub struct ClusterOptions {
    pub cutoff1: f64,
    pub cutoff2: f64,
    pub intron_weight: f64,
    pub sw_score: i64,
    pub name2_mode: Name2Mode,
    pub batch_size: usize,
    pub batch_rounds: usize,
}

impl Default for ClusterOptions {
    fn default() -> Self {
        Self {
            cutoff1: DEFAULT_CUTOFF1,
            cutoff2: DEFAULT_CUTOFF2,
            intron_weight: DEFAULT_INTRON_WEIGHT,
            sw_score: DEFAULT_SW_SCORE,
            name2_mode: Name2Mode::Full,
            batch_size: 0,
            batch_rounds: 100,
        }
    }
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

fn merge_tracks_by_name(tracks: Vec<Track>) -> Vec<Track> {
    let mut merged: Vec<Track> = Vec::with_capacity(tracks.len());
    let mut positions: HashMap<String, usize> = HashMap::new();

    for track in tracks {
        if let Some(&idx) = positions.get(track.tx.name.as_str()) {
            merged[idx].subreads.extend(track.subreads);
            continue;
        }

        positions.insert(track.tx.name.clone(), merged.len());
        merged.push(track);
    }

    merged
}

fn split_reference_and_read_tracks(tracks: Vec<Track>) -> (Vec<Track>, Vec<Track>) {
    let mut refs: Vec<Track> = Vec::new();
    let mut reads: Vec<Track> = Vec::new();

    for track in tracks {
        if is_isoform_anno(&track.tx) {
            refs.push(track);
        } else {
            reads.push(track);
        }
    }

    (refs, reads)
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
    sw_score: i64,
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
                    if mode == DistanceMode::Ratio || i64::from(tracks[i].tx.score) < sw_score {
                        drop.insert(i);
                        merge_subreads(i, j, &mut tracks);
                    }
                }
                std::cmp::Ordering::Equal => {
                    drop.insert(i);
                    merge_subreads(i, j, &mut tracks);
                }
                std::cmp::Ordering::Greater => {
                    if mode == DistanceMode::Ratio || i64::from(tracks[j].tx.score) < sw_score {
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

fn update_name2(isoforms: &mut [Track], ref_names: &HashSet<String>, mode: Name2Mode) {
    if mode == Name2Mode::None {
        for track in isoforms.iter_mut() {
            set_extra(&mut track.tx, NAME2_COL, "none".to_owned());
        }
        return;
    }

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
                let mut coverage = 0.0f64;
                for name in &track.subreads {
                    if ref_names.contains(name) {
                        continue;
                    }
                    let denom = occurrence.get(name.as_str()).copied().unwrap_or(0);
                    if denom > 0 {
                        coverage += 1.0f64 / denom as f64;
                    }
                }

                match mode {
                    Name2Mode::Full => {
                        let mut subreads: Vec<&str> =
                            track.subreads.iter().map(|s| s.as_str()).collect();
                        subreads.sort_unstable();
                        let joined = subreads.join(",");
                        format!("{joined},|{coverage}")
                    }
                    Name2Mode::Coverage => format!("|{coverage}"),
                    Name2Mode::None => unreachable!("handled above"),
                }
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

fn cluster_once(tracks: Vec<Track>, options: ClusterOptions) -> Vec<Track> {
    let tracks = merge_tracks_by_name(tracks);
    let tracks = filter_pass(
        tracks,
        DistanceMode::Ratio,
        options.cutoff1,
        options.intron_weight,
        options.sw_score,
    );
    filter_pass(
        tracks,
        DistanceMode::RatioShort,
        options.cutoff2,
        options.intron_weight,
        options.sw_score,
    )
}

fn batch_overlap_merge(tracks: Vec<Track>, options: ClusterOptions) -> Vec<Track> {
    if options.batch_size == 0 {
        return cluster_once(tracks, options);
    }

    let batch_size = options.batch_size.max(1);
    let max_rounds = options.batch_rounds;

    let tracks = merge_tracks_by_name(tracks);
    let (mut anchors, mut pending_reads) = split_reference_and_read_tracks(tracks);

    let mut rounds = 0usize;
    while rounds < max_rounds && pending_reads.len() > batch_size {
        let remainder = pending_reads.split_off(batch_size);
        let batch = std::mem::replace(&mut pending_reads, remainder);

        let mut combined: Vec<Track> = Vec::with_capacity(anchors.len() + batch.len());
        combined.extend(anchors.iter().cloned());
        combined.extend(batch);

        let merged_batch = cluster_once(combined, options);
        let (next_anchors, batch_reads) = split_reference_and_read_tracks(merged_batch);

        anchors = merge_tracks_by_name(next_anchors);

        let mut next_pending: Vec<Track> =
            Vec::with_capacity(batch_reads.len() + pending_reads.len());
        next_pending.extend(batch_reads);
        next_pending.extend(std::mem::take(&mut pending_reads));
        pending_reads = merge_tracks_by_name(next_pending);

        rounds += 1;
    }

    let mut combined: Vec<Track> = Vec::with_capacity(anchors.len() + pending_reads.len());
    combined.extend(anchors);
    combined.extend(pending_reads);
    cluster_once(combined, options)
}

fn process_partition(
    references: &[Transcript],
    reads: &[Transcript],
    ref_indices: &[usize],
    read_indices: &[usize],
    ref_names: &HashSet<String>,
    options: ClusterOptions,
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

        let mut tracks = batch_overlap_merge(tracks, options);

        update_name2(&mut tracks, ref_names, options.name2_mode);
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
    cluster_with_options(reads, references, threads, ClusterOptions::default())
}

pub fn cluster_with_options(
    reads: &[Transcript],
    references: Option<&[Transcript]>,
    threads: usize,
    options: ClusterOptions,
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
                options,
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
                        options,
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

    #[test]
    fn name2_mode_coverage_writes_only_coverage_but_keeps_mapping() {
        let refs = vec![make_tx(
            "ref",
            Strand::Plus,
            &[(100, 110), (120, 130), (140, 150)],
            "isoform_anno",
            100,
        )];
        let reads = vec![make_tx(
            "read1",
            Strand::Plus,
            &[(100, 110), (120, 130), (140, 150)],
            "nanopore_read",
            0,
        )];

        let result = cluster_with_options(
            &reads,
            Some(&refs),
            1,
            ClusterOptions {
                name2_mode: Name2Mode::Coverage,
                ..ClusterOptions::default()
            },
        );

        assert_eq!(result.read_to_isoform.len(), 1);
        let ref_isoform = result
            .isoforms
            .iter()
            .find(|tx| tx.name == "ref")
            .expect("reference isoform retained");
        let name2 = ref_isoform
            .extra_fields
            .first()
            .expect("name2 payload present");
        assert!(name2.starts_with('|'));
        assert!(!name2.contains("read1"));
    }

    #[test]
    fn sw_score_minus_one_disables_ratio_short_truncation_merge() {
        use std::collections::HashSet;

        let refs = vec![make_tx(
            "ref",
            Strand::Plus,
            &[(100, 110), (120, 130), (140, 150)],
            "isoform_anno",
            100,
        )];
        let reads = vec![
            make_tx(
                "read_long",
                Strand::Plus,
                &[(100, 110), (120, 130), (140, 150)],
                "nanopore_read",
                20,
            ),
            make_tx(
                "read_short",
                Strand::Plus,
                &[(120, 130), (140, 150)],
                "nanopore_read",
                0,
            ),
        ];

        let merged = cluster_with_options(&reads, Some(&refs), 1, ClusterOptions::default());
        let no_merge = cluster_with_options(
            &reads,
            Some(&refs),
            1,
            ClusterOptions {
                sw_score: -1,
                ..ClusterOptions::default()
            },
        );

        let merged_reads: HashSet<&str> = merged
            .read_to_isoform
            .iter()
            .map(|(read_id, _)| read_id.as_str())
            .collect();
        let no_merge_reads: HashSet<&str> = no_merge
            .read_to_isoform
            .iter()
            .map(|(read_id, _)| read_id.as_str())
            .collect();
        assert_eq!(merged_reads.len(), 2);
        assert_eq!(no_merge_reads.len(), 2);

        let no_merge_targets: HashSet<&str> = no_merge
            .read_to_isoform
            .iter()
            .map(|(_, isoform_id)| isoform_id.as_str())
            .collect();

        assert_eq!(no_merge_targets.len(), 2);

        let merged_short_targets: HashSet<&str> = merged
            .read_to_isoform
            .iter()
            .filter(|(read_id, _)| read_id == "read_short")
            .map(|(_, isoform_id)| isoform_id.as_str())
            .collect();
        let no_merge_short_targets: HashSet<&str> = no_merge
            .read_to_isoform
            .iter()
            .filter(|(read_id, _)| read_id == "read_short")
            .map(|(_, isoform_id)| isoform_id.as_str())
            .collect();

        assert_eq!(merged_short_targets.len(), 2);
        assert!(merged_short_targets.contains("read_long"));
        assert!(merged_short_targets.contains("ref"));
        assert_eq!(no_merge_short_targets.len(), 1);
        assert!(no_merge_short_targets.contains("read_short"));
    }

    #[test]
    fn sw_score_equal_to_cutoff_keeps_sl_read_as_its_own_track() {
        use std::collections::HashSet;

        let refs = vec![make_tx(
            "ref",
            Strand::Plus,
            &[(100, 110), (120, 130), (140, 150)],
            "isoform_anno",
            100,
        )];
        let reads = vec![
            make_tx(
                "read_long",
                Strand::Plus,
                &[(100, 110), (120, 130), (140, 150)],
                "nanopore_read",
                20,
            ),
            make_tx(
                "read_sl",
                Strand::Plus,
                &[(120, 130), (140, 150)],
                "nanopore_read",
                DEFAULT_SW_SCORE as u32,
            ),
        ];

        let result = cluster_with_options(&reads, Some(&refs), 1, ClusterOptions::default());

        let sl_targets: HashSet<&str> = result
            .read_to_isoform
            .iter()
            .filter(|(read_id, _)| read_id == "read_sl")
            .map(|(_, isoform_id)| isoform_id.as_str())
            .collect();

        assert_eq!(sl_targets.len(), 1);
        assert!(sl_targets.contains("read_sl"));
    }

    #[test]
    fn overlap_batching_matches_unbatched_on_simple_locus() {
        let refs = vec![make_tx(
            "ref",
            Strand::Plus,
            &[(100, 110), (120, 130), (140, 150)],
            "isoform_anno",
            100,
        )];
        let reads = vec![
            make_tx(
                "read1",
                Strand::Plus,
                &[(100, 110), (120, 130), (140, 150)],
                "nanopore_read",
                0,
            ),
            make_tx(
                "read2",
                Strand::Plus,
                &[(100, 110), (120, 130), (140, 150)],
                "nanopore_read",
                0,
            ),
            make_tx(
                "read3",
                Strand::Plus,
                &[(120, 130), (140, 150)],
                "nanopore_read",
                0,
            ),
        ];

        let single_pass = cluster_with_options(&reads, Some(&refs), 1, ClusterOptions::default());
        let batched = cluster_with_options(
            &reads,
            Some(&refs),
            1,
            ClusterOptions {
                batch_size: 1,
                batch_rounds: 10,
                ..ClusterOptions::default()
            },
        );

        assert_eq!(single_pass.read_to_isoform, batched.read_to_isoform);

        let single_names: Vec<&str> = single_pass
            .isoforms
            .iter()
            .map(|tx| tx.name.as_str())
            .collect();
        let batch_names: Vec<&str> = batched.isoforms.iter().map(|tx| tx.name.as_str()).collect();
        assert_eq!(single_names, batch_names);
    }
}
