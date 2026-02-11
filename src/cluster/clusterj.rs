use std::collections::{HashMap, HashSet};

use crate::cluster::result::ClusterResult;
use crate::model::{Coord, Interval, Strand, Transcript};

#[derive(Clone, Debug)]
struct Track {
    tx: Transcript,
    subreads: HashSet<String>,
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct PartitionKey {
    chrom: String,
    strand: Strand,
}

const NAME2_COL: usize = 0;
const TTYPE_COL: usize = 4;
const DEFAULT_SW_SCORE: i64 = 11;

fn get_extra(tx: &Transcript, idx: usize) -> Option<&str> {
    tx.extra_fields.get(idx).map(|value| value.as_str())
}

fn set_extra(tx: &mut Transcript, idx: usize, value: String) {
    if tx.extra_fields.len() <= idx {
        tx.extra_fields.resize(idx + 1, "none".to_owned());
    }
    tx.extra_fields[idx] = value;
}

fn ttype(tx: &Transcript) -> Option<&str> {
    get_extra(tx, TTYPE_COL)
}

fn is_isoform_anno(tx: &Transcript) -> bool {
    matches!(ttype(tx), Some("isoform_anno"))
}

fn track_weight(tx: &Transcript, ref_weight: u32, read_weight: u32) -> u32 {
    match ttype(tx) {
        Some("isoform_anno") => ref_weight,
        Some("nanopore_read") => read_weight,
        _ => 1,
    }
}

fn exon_len(tx: &Transcript) -> u32 {
    tx.exons.iter().map(|exon| exon.len()).sum()
}

fn junction_positions(tx: &Transcript) -> Vec<u32> {
    let mut boundaries: Vec<u32> = Vec::with_capacity(tx.exons.len() * 2);
    for exon in &tx.exons {
        boundaries.push(exon.start.get());
        boundaries.push(exon.end.get());
    }

    match tx.strand {
        Strand::Plus | Strand::Unknown => {}
        Strand::Minus => boundaries.reverse(),
    }

    if boundaries.len() <= 2 {
        return Vec::new();
    }
    boundaries[1..boundaries.len() - 1].to_vec()
}

fn rebuild_exons_from_junctions(
    tx_start: Coord,
    tx_end: Coord,
    junctions: &[u32],
) -> Vec<Interval> {
    let start = tx_start.get();
    let end = tx_end.get();

    let mut junctions = junctions.to_vec();
    junctions.sort_unstable();

    // Junction correction can drift slightly outside the transcript span; clamp defensively so we
    // never panic when rebuilding exon intervals.
    for pos in &mut junctions {
        *pos = (*pos).clamp(start, end);
    }
    junctions.dedup();

    if junctions.is_empty() || !junctions.len().is_multiple_of(2) || start >= end {
        return vec![Interval::new(tx_start, tx_end).expect("valid span")];
    }

    let mut exons: Vec<Interval> = Vec::new();
    exons.push(Interval::new(Coord::new(start), Coord::new(junctions[0])).expect("valid exon"));

    let mut idx = 1usize;
    while idx + 1 < junctions.len() {
        let exon = Interval::new(Coord::new(junctions[idx]), Coord::new(junctions[idx + 1]))
            .expect("valid exon");
        exons.push(exon);
        idx += 2;
    }

    exons.push(
        Interval::new(
            Coord::new(*junctions.last().expect("non-empty")),
            Coord::new(end),
        )
        .expect("valid exon"),
    );

    exons
}

fn group_consecutive_indices(indices: &[usize]) -> Vec<Vec<usize>> {
    let mut groups: Vec<Vec<usize>> = Vec::new();
    for &idx in indices {
        match groups.last_mut() {
            Some(group) if idx == group.last().copied().unwrap_or(idx).saturating_add(1) => {
                group.push(idx);
            }
            _ => groups.push(vec![idx]),
        }
    }
    groups
}

fn build_corrected_site_maps(
    site_cov: &HashMap<u32, u32>,
    cov_cutoff: u32,
    pos_cutoff: u32,
) -> (HashMap<u32, u32>, HashSet<u32>) {
    let mut keys: Vec<u32> = site_cov.keys().copied().collect();
    keys.sort_unstable();

    let mut high_by_idx: HashMap<usize, u32> = HashMap::new();
    let mut low_indices: Vec<usize> = Vec::new();
    let mut low_by_idx: HashMap<usize, u32> = HashMap::new();

    for (idx, coord) in keys.iter().enumerate() {
        let cov = site_cov.get(coord).copied().unwrap_or(0);
        if cov >= cov_cutoff {
            high_by_idx.insert(idx, *coord);
        } else {
            low_indices.push(idx);
            low_by_idx.insert(idx, *coord);
        }
    }

    let low_groups = group_consecutive_indices(&low_indices);
    let mut w_to_r: HashMap<u32, u32> = HashMap::new();
    let mut w_to_no: HashSet<u32> = HashSet::new();

    for group in low_groups {
        let start = *group.first().expect("non-empty group");
        let end = *group.last().expect("non-empty group");
        for site_idx in group {
            let low_coord = low_by_idx.get(&site_idx).copied().unwrap_or(0);

            let corrected = high_by_idx
                .get(&(end + 1))
                .copied()
                .filter(|high_coord| low_coord.abs_diff(*high_coord) <= pos_cutoff)
                .or_else(|| {
                    if start == 0 {
                        return None;
                    }
                    high_by_idx
                        .get(&(start - 1))
                        .copied()
                        .filter(|high_coord| low_coord.abs_diff(*high_coord) <= pos_cutoff)
                });

            if let Some(corrected_coord) = corrected {
                w_to_r.insert(low_coord, corrected_coord);
            } else {
                w_to_no.insert(low_coord);
            }
        }
    }

    (w_to_r, w_to_no)
}

fn flow_junction_correct(
    tracks: Vec<Track>,
    cov_cutoff: u32,
    pos_cutoff: u32,
) -> (Vec<Track>, Vec<Track>) {
    let junctions_cache: Vec<Vec<u32>> = tracks
        .iter()
        .map(|track| junction_positions(&track.tx))
        .collect();

    let mut site_cov: HashMap<u32, u32> = HashMap::new();
    for (track, junctions) in tracks.iter().zip(junctions_cache.iter()) {
        let weight = track_weight(&track.tx, 5, 1);
        for &pos in junctions {
            *site_cov.entry(pos).or_insert(0) += weight;
        }
    }

    let (w_to_r, w_to_no) = build_corrected_site_maps(&site_cov, cov_cutoff, pos_cutoff);

    let mut corrected: Vec<Track> = Vec::new();
    let mut rare: Vec<Track> = Vec::new();

    for (idx, mut track) in tracks.into_iter().enumerate() {
        let junctions = &junctions_cache[idx];
        if junctions.iter().any(|pos| w_to_no.contains(pos)) {
            rare.push(track);
            continue;
        }

        if !is_isoform_anno(&track.tx) && !junctions.is_empty() {
            let mut corrected_junctions: Vec<u32> = Vec::with_capacity(junctions.len());
            let mut changed = false;
            for &pos in junctions {
                let corrected_pos = w_to_r.get(&pos).copied().unwrap_or(pos);
                if corrected_pos != pos {
                    changed = true;
                }
                corrected_junctions.push(corrected_pos);
            }

            if changed && track.tx.extra_fields.len() > TTYPE_COL {
                set_extra(
                    &mut track.tx,
                    TTYPE_COL,
                    "nanopore_read_corrected".to_owned(),
                );
            }

            track.tx.exons = rebuild_exons_from_junctions(
                track.tx.tx_start,
                track.tx.tx_end,
                &corrected_junctions,
            );
        }

        corrected.push(track);
    }

    (corrected, rare)
}

fn junctions_equal(a: &[u32], b: &[u32], offset: u32) -> bool {
    if offset == 0 {
        return a == b;
    }

    let mut matched_a: HashSet<u32> = HashSet::new();
    let mut matched_b: HashSet<u32> = HashSet::new();

    for &i in a {
        for &j in b {
            if i.abs_diff(j) <= offset {
                matched_a.insert(i);
                matched_b.insert(j);
            }
        }
    }

    a.iter().copied().collect::<HashSet<u32>>() == matched_a
        && b.iter().copied().collect::<HashSet<u32>>() == matched_b
}

fn fuzzy_intersection(a: &[u32], b: &[u32], offset: u32) -> HashMap<u32, u32> {
    let mut match_dic: HashMap<u32, u32> = HashMap::new();
    for &i in a {
        for &j in b {
            if i.abs_diff(j) <= offset {
                match_dic.insert(i, j);
            }
        }
    }
    match_dic
}

fn compare_ei_by_boundary(a: &[u32], reference: &[u32], offset: u32) -> (Vec<usize>, Vec<usize>) {
    if offset == 0 {
        let ascending = match (reference.first(), reference.last()) {
            (Some(first), Some(last)) => first <= last,
            _ => true,
        };

        let mut missed_order: Vec<usize> = Vec::new();
        let mut extra_order: Vec<usize> = Vec::new();

        let mut i = 0usize;
        let mut j = 0usize;
        while i < a.len() && j < reference.len() {
            let ai = a[i];
            let rj = reference[j];
            if ai == rj {
                i += 1;
                j += 1;
                continue;
            }

            let a_before_ref = if ascending { ai < rj } else { ai > rj };
            if a_before_ref {
                extra_order.push(i);
                i += 1;
            } else {
                missed_order.push(j);
                j += 1;
            }
        }

        while i < a.len() {
            extra_order.push(i);
            i += 1;
        }
        while j < reference.len() {
            missed_order.push(j);
            j += 1;
        }

        return (missed_order, extra_order);
    }

    let match_dic = fuzzy_intersection(a, reference, offset);
    let junction_new: Vec<u32> = a
        .iter()
        .copied()
        .map(|pos| match_dic.get(&pos).copied().unwrap_or(pos))
        .collect();

    let mut posdic_a: HashMap<u32, usize> = HashMap::new();
    for (idx, pos) in junction_new.iter().copied().enumerate() {
        posdic_a.insert(pos, idx);
    }

    let mut posdic_ref: HashMap<u32, usize> = HashMap::new();
    for (idx, pos) in reference.iter().copied().enumerate() {
        posdic_ref.insert(pos, idx);
    }

    let set_a: HashSet<u32> = junction_new.iter().copied().collect();
    let set_ref: HashSet<u32> = reference.iter().copied().collect();

    let mut missed_order: Vec<usize> = set_ref
        .difference(&set_a)
        .filter_map(|pos| posdic_ref.get(pos).copied())
        .collect();
    missed_order.sort_unstable();

    let mut extra_order: Vec<usize> = set_a
        .difference(&set_ref)
        .filter_map(|pos| posdic_a.get(pos).copied())
        .collect();
    extra_order.sort_unstable();

    (missed_order, extra_order)
}

fn is_junction_5primer(missed_order: &[usize]) -> bool {
    if missed_order.is_empty() || missed_order[0] != 0 {
        return false;
    }

    let groups = group_consecutive_indices(missed_order);
    groups.len() == 1
}

fn is_junction_inside(
    short_junctions: &[u32],
    short_exon_len: u32,
    long_junctions: &[u32],
    long_exon_len: u32,
) -> bool {
    if short_junctions.is_empty() {
        return false;
    }

    if junctions_equal(short_junctions, long_junctions, 0) {
        return short_exon_len < long_exon_len;
    }

    let (missed_order, extra_order) = compare_ei_by_boundary(short_junctions, long_junctions, 0);
    if missed_order.is_empty() {
        return false;
    }
    is_junction_5primer(&missed_order) && extra_order.is_empty()
}

fn is_single_exon_in(single: &Transcript, other: &Transcript, other_junctions: &[u32]) -> bool {
    if other_junctions.is_empty() {
        match single.strand {
            Strand::Plus => {
                if single.tx_start.get() <= other.tx_start.get()
                    || single.tx_end.get() >= other.tx_end.get()
                {
                    return false;
                }
                true
            }
            Strand::Minus => {
                if single.tx_start.get() >= other.tx_start.get()
                    || single.tx_end.get() <= other.tx_end.get()
                {
                    return false;
                }
                true
            }
            Strand::Unknown => false,
        }
    } else {
        let last_junction = *other_junctions.last().expect("non-empty");
        match single.strand {
            Strand::Plus => {
                if single.tx_start.get() <= last_junction
                    || single.tx_end.get() >= other.tx_end.get()
                {
                    return false;
                }
                true
            }
            Strand::Minus => {
                if single.tx_start.get() <= other.tx_start.get()
                    || single.tx_end.get() >= last_junction
                {
                    return false;
                }
                true
            }
            Strand::Unknown => false,
        }
    }
}

fn get_two_mut<T>(slice: &mut [T], i: usize, j: usize) -> (&mut T, &mut T) {
    assert!(i != j, "indices must be distinct");
    if i < j {
        let (left, right) = slice.split_at_mut(j);
        (&mut left[i], &mut right[0])
    } else {
        let (left, right) = slice.split_at_mut(i);
        (&mut right[0], &mut left[j])
    }
}

fn junction_simple_merge(tracks: &mut [Track], sw_score: i64) -> Vec<usize> {
    let junctions_cache: Vec<Vec<u32>> = tracks
        .iter()
        .map(|track| junction_positions(&track.tx))
        .collect();
    let exon_lens: Vec<u32> = tracks.iter().map(|track| exon_len(&track.tx)).collect();

    let mut dropped: Vec<bool> = vec![false; tracks.len()];
    for i in 0..tracks.len() {
        if dropped[i] {
            continue;
        }

        let score_i = i64::from(tracks[i].tx.score);
        if score_i > sw_score {
            continue;
        }

        let junctions_i = &junctions_cache[i];
        let exon_len_i = exon_lens[i];
        for j in 0..tracks.len() {
            if i == j {
                continue;
            }
            if dropped[j] && !is_isoform_anno(&tracks[j].tx) {
                continue;
            }

            if junctions_i.is_empty() {
                if is_single_exon_in(&tracks[i].tx, &tracks[j].tx, &junctions_cache[j]) {
                    dropped[i] = true;
                    let (short, long) = get_two_mut(tracks, i, j);
                    let name = short.tx.name.clone();
                    long.subreads.insert(name);
                    long.subreads.extend(short.subreads.iter().cloned());
                }
            } else if is_junction_inside(junctions_i, exon_len_i, &junctions_cache[j], exon_lens[j])
            {
                dropped[i] = true;
                let (short, long) = get_two_mut(tracks, i, j);
                let name = short.tx.name.clone();
                long.subreads.insert(name);
                long.subreads.extend(short.subreads.iter().cloned());
            }
        }
    }

    let mut keep_vec: Vec<usize> = Vec::with_capacity(tracks.len());
    for (idx, track) in tracks.iter().enumerate() {
        if !dropped[idx] || is_isoform_anno(&track.tx) {
            keep_vec.push(idx);
        }
    }
    keep_vec
}

fn select_tracks_by_keep_indices(tracks: Vec<Track>, keep_indices: Vec<usize>) -> Vec<Track> {
    if keep_indices.len() == tracks.len() {
        return tracks;
    }

    let mut keep_mask = vec![false; tracks.len()];
    for idx in keep_indices {
        if idx < keep_mask.len() {
            keep_mask[idx] = true;
        }
    }

    tracks
        .into_iter()
        .enumerate()
        .filter_map(|(idx, track)| keep_mask[idx].then_some(track))
        .collect()
}

fn merge_tracks_by_name(tracks: Vec<Track>) -> Vec<Track> {
    let mut out: Vec<Track> = Vec::new();
    let mut index_by_name: HashMap<String, usize> = HashMap::new();

    for track in tracks {
        if let Some(&idx) = index_by_name.get(&track.tx.name) {
            out[idx].subreads.extend(track.subreads);
        } else {
            index_by_name.insert(track.tx.name.clone(), out.len());
            out.push(track);
        }
    }

    out
}

fn batch_junction_simple_merge(
    tracks: Vec<Track>,
    sw_score: i64,
    batch_size: usize,
    max_rounds: usize,
) -> Vec<Track> {
    let batch_size = batch_size.max(1);
    let max_rounds = max_rounds.max(1);

    let mut tracks = tracks;
    let mut rounds = 0usize;
    let mut previous_len = tracks.len();

    // Match TrackCluster Python `process_one_junction_corrected_try` behavior:
    // - Intermediate merges always use the default SW cutoff (11) even if `sw_score` is -1.
    // - The final merge uses the user-provided `sw_score`.
    while rounds < max_rounds && tracks.len() > batch_size {
        let rest = tracks.split_off(batch_size);
        let mut batch = tracks;

        let keep_indices = junction_simple_merge(&mut batch, DEFAULT_SW_SCORE);
        let kept_batch = select_tracks_by_keep_indices(batch, keep_indices);

        let mut combined: Vec<Track> = Vec::with_capacity(kept_batch.len() + rest.len());
        combined.extend(kept_batch);
        combined.extend(rest);
        tracks = merge_tracks_by_name(combined);

        if tracks.len() >= previous_len {
            break;
        }
        previous_len = tracks.len();
        rounds += 1;
    }

    let keep_indices = junction_simple_merge(&mut tracks, sw_score);
    select_tracks_by_keep_indices(tracks, keep_indices)
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
    unused: Vec<Transcript>,
}

struct WorkItem {
    index: usize,
    refs: Vec<Transcript>,
    reads: Vec<Transcript>,
}

fn process_partition(
    refs: Vec<Transcript>,
    reads: Vec<Transcript>,
    ref_names: &HashSet<String>,
    sw_score: i64,
    batch_size: usize,
    batch_rounds: usize,
) -> PartitionResult {
    let mut tracks: Vec<Track> = Vec::with_capacity(refs.len() + reads.len());
    for tx in refs.into_iter().chain(reads.into_iter()) {
        tracks.push(Track {
            tx,
            subreads: HashSet::new(),
        });
    }

    let (mut corrected, rare) = flow_junction_correct(tracks, 2, 10);
    let unused: Vec<Transcript> = rare.into_iter().map(|track| track.tx).collect();

    let mut kept = if batch_size == 0 {
        let keep_indices = junction_simple_merge(&mut corrected, sw_score);
        select_tracks_by_keep_indices(corrected, keep_indices)
    } else {
        batch_junction_simple_merge(corrected, sw_score, batch_size, batch_rounds)
    };

    update_name2_with_coverage(&mut kept, ref_names);
    let pairs = build_read_to_isoform(&kept, ref_names);
    let isoforms = kept.into_iter().map(|track| track.tx).collect();

    PartitionResult {
        isoforms,
        pairs,
        unused,
    }
}

pub fn clusterj(
    reads: &[Transcript],
    references: Option<&[Transcript]>,
    threads: usize,
    sw_score: i64,
    batch_size: usize,
    batch_rounds: usize,
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

    let mut refs_by_key: HashMap<PartitionKey, Vec<Transcript>> = HashMap::new();
    for tx in references {
        refs_by_key
            .entry(PartitionKey {
                chrom: tx.chrom.clone(),
                strand: tx.strand,
            })
            .or_default()
            .push(tx.clone());
    }

    let mut reads_by_key: HashMap<PartitionKey, Vec<Transcript>> = HashMap::new();
    for read in reads {
        let key = PartitionKey {
            chrom: read.chrom.clone(),
            strand: read.strand,
        };
        if !refs_by_key.contains_key(&key) {
            continue;
        }
        reads_by_key.entry(key).or_default().push(read.clone());
    }

    let mut all_isoforms: Vec<Transcript> = Vec::new();
    let mut all_pairs: Vec<(String, String)> = Vec::new();
    let mut all_unused: Vec<Transcript> = Vec::new();

    let mut keys: Vec<PartitionKey> = refs_by_key.keys().cloned().collect();
    keys.sort_by(|a, b| a.chrom.cmp(&b.chrom).then_with(|| a.strand.cmp(&b.strand)));

    let mut work: Vec<WorkItem> = Vec::with_capacity(keys.len());
    for (index, key) in keys.iter().enumerate() {
        work.push(WorkItem {
            index,
            refs: refs_by_key.remove(key).unwrap_or_default(),
            reads: reads_by_key.remove(key).unwrap_or_default(),
        });
    }

    let mut parts: Vec<Option<PartitionResult>> = (0..keys.len()).map(|_| None).collect();
    if threads == 1 || work.len() <= 1 {
        for item in work {
            parts[item.index] = Some(process_partition(
                item.refs,
                item.reads,
                &ref_names,
                sw_score,
                batch_size,
                batch_rounds,
            ));
        }
    } else {
        use std::sync::{mpsc, Arc, Mutex};

        let queue = Arc::new(Mutex::new(work));
        let (tx, rx) = mpsc::channel::<(usize, PartitionResult)>();

        let worker_count = threads.min(keys.len());
        let mut handles = Vec::with_capacity(worker_count);
        for _ in 0..worker_count {
            let queue = Arc::clone(&queue);
            let tx = tx.clone();
            let ref_names = Arc::clone(&ref_names);

            handles.push(std::thread::spawn(move || loop {
                let item = {
                    let mut guard = queue.lock().expect("work queue poisoned");
                    guard.pop()
                };
                let Some(item) = item else {
                    break;
                };

                let result = process_partition(
                    item.refs,
                    item.reads,
                    &ref_names,
                    sw_score,
                    batch_size,
                    batch_rounds,
                );
                if tx.send((item.index, result)).is_err() {
                    break;
                }
            }));
        }
        drop(tx);

        for _ in 0..keys.len() {
            let (index, result) = rx.recv().expect("worker hung up");
            parts[index] = Some(result);
        }

        for handle in handles {
            handle.join().expect("worker panicked");
        }
    }

    for part in parts.into_iter().flatten() {
        all_isoforms.extend(part.isoforms);
        all_pairs.extend(part.pairs);
        all_unused.extend(part.unused);
    }

    ClusterResult {
        isoforms: all_isoforms,
        read_to_isoform: all_pairs,
        unused: all_unused,
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
                    name.to_owned(),
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
    fn merges_5prime_truncations_into_multiple_compatible_tracks() {
        let refs = vec![
            make_tx(
                "ref_a",
                Strand::Plus,
                &[(100, 110), (120, 130), (140, 150)],
                "isoform_anno",
                100,
            ),
            make_tx(
                "ref_b",
                Strand::Plus,
                &[(80, 90), (100, 110), (120, 130), (140, 150)],
                "isoform_anno",
                100,
            ),
        ];
        let reads = vec![make_tx(
            "read_trunc",
            Strand::Plus,
            &[(120, 130), (140, 150)],
            "nanopore_read",
            1,
        )];

        let result = clusterj(&reads, Some(&refs), 1, 11, 0, 1);
        assert_eq!(result.unused.len(), 0);

        let iso_names: HashSet<_> = result.isoforms.iter().map(|tx| tx.name.as_str()).collect();
        assert!(iso_names.contains("ref_a"));
        assert!(iso_names.contains("ref_b"));
        assert!(!iso_names.contains("read_trunc"));

        let subread_sets: HashMap<&str, HashSet<&str>> = result
            .isoforms
            .iter()
            .map(|tx| {
                let raw = tx.extra_fields.first().map(|s| s.as_str()).unwrap_or("");
                let sub = raw.split(",|").next().unwrap_or("");
                let set: HashSet<&str> = sub.split(',').filter(|s| !s.is_empty()).collect();
                (tx.name.as_str(), set)
            })
            .collect();

        assert!(subread_sets.get("ref_a").unwrap().contains("read_trunc"));
        assert!(subread_sets.get("ref_b").unwrap().contains("read_trunc"));
    }

    #[test]
    fn threaded_execution_is_deterministic() {
        let refs = vec![
            make_tx(
                "ref_plus",
                Strand::Plus,
                &[(0, 10), (20, 30)],
                "isoform_anno",
                100,
            ),
            make_tx(
                "ref_minus",
                Strand::Minus,
                &[(100, 110), (120, 130)],
                "isoform_anno",
                100,
            ),
        ];
        let reads = vec![
            make_tx("read_plus", Strand::Plus, &[(20, 30)], "nanopore_read", 1),
            make_tx(
                "read_minus",
                Strand::Minus,
                &[(120, 130)],
                "nanopore_read",
                1,
            ),
        ];

        let single = clusterj(&reads, Some(&refs), 1, 11, 0, 1);
        let threaded = clusterj(&reads, Some(&refs), 4, 11, 0, 1);
        assert_eq!(single.isoforms, threaded.isoforms);
        assert_eq!(single.read_to_isoform, threaded.read_to_isoform);
        assert_eq!(single.unused, threaded.unused);
    }
}
