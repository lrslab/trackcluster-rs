use std::collections::{HashMap, HashSet};

use crate::cluster::{clusterj::Name2Mode, result::ClusterResult};
use crate::interval::{cluster_by_span, exonic_overlap_bp, StrandMode};
use crate::model::{Interval, Strand, Transcript};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
enum TrackSource {
    Reference,
    Read,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
struct TrackPair {
    left: usize,
    right: usize,
}

impl TrackPair {
    fn new(a: usize, b: usize) -> Self {
        if a <= b {
            Self { left: a, right: b }
        } else {
            Self { left: b, right: a }
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct PairOverlap {
    pair: TrackPair,
    exon_overlap: u32,
}

#[derive(Clone, Copy, Debug)]
struct ExonRecord {
    track_idx: usize,
    interval: Interval,
}

#[derive(Clone, Copy, Debug)]
struct FilterParams {
    mode: DistanceMode,
    cutoff: f64,
    intron_weight: f64,
    sw_score: i64,
}

#[derive(Clone, Copy, Debug)]
struct PairDistance {
    left: usize,
    right: usize,
    exon_overlap: u32,
    intron_overlap: u32,
}

#[derive(Clone, Debug)]
struct Track {
    tx: Transcript,
    source: TrackSource,
    subreads: HashSet<ReadInstance>,
    exon_len: u32,
    introns: Vec<Interval>,
    intron_len: u32,
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct ReadInstance {
    index: usize,
    name: String,
}

impl Track {
    fn new(tx: Transcript, source: TrackSource, read_index: Option<usize>) -> Self {
        let mut subreads = HashSet::new();
        if source == TrackSource::Read {
            subreads.insert(ReadInstance {
                index: read_index.expect("read tracks require an input-instance index"),
                name: tx.name.clone(),
            });
        } else {
            debug_assert!(read_index.is_none());
        }
        let exon_len = tx.exons.iter().map(|exon| exon.len()).sum();
        let introns = tx.introns();
        let intron_len = introns.iter().map(|intron| intron.len()).sum();
        Self {
            tx,
            source,
            subreads,
            exon_len,
            introns,
            intron_len,
        }
    }

    fn reference(tx: Transcript) -> Self {
        Self::new(tx, TrackSource::Reference, None)
    }

    fn read(tx: Transcript, index: usize) -> Self {
        Self::new(tx, TrackSource::Read, Some(index))
    }

    fn is_reference(&self) -> bool {
        self.source == TrackSource::Reference
    }

    fn is_read(&self) -> bool {
        self.source == TrackSource::Read
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
const SPARSE_PAIR_THRESHOLD: usize = 512;

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

impl ClusterOptions {
    /// Validate scientific overlap options.
    pub fn validate(self) -> Result<(), crate::config::ParameterError> {
        crate::config::BatchRounds::new(self.batch_rounds)?;
        crate::config::UnitFraction::new("overlap pass-1 cutoff", self.cutoff1)?;
        crate::config::UnitFraction::new("overlap pass-2 cutoff", self.cutoff2)?;
        crate::config::NonNegativeWeight::new("overlap intron weight", self.intron_weight)?;
        Ok(())
    }
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

fn combined_distance_from_overlaps(
    a: &Track,
    b: &Track,
    mode: DistanceMode,
    intron_weight: f64,
    exon_overlap: u32,
    intron_overlap: u32,
) -> f64 {
    let exon_overlap = exon_overlap as f64;
    let a_exon_len = a.exon_len as f64;
    let b_exon_len = b.exon_len as f64;

    let exon_denom = match mode {
        DistanceMode::Ratio => (a_exon_len + b_exon_len - exon_overlap).max(1.0),
        DistanceMode::RatioShort => a_exon_len.min(b_exon_len).max(1.0),
    };
    let exon_dist = 1.0 - (exon_overlap / exon_denom);

    let intron_overlap = intron_overlap as f64;
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

fn exon_record_cmp(left: &ExonRecord, right: &ExonRecord, tracks: &[Track]) -> std::cmp::Ordering {
    let left_tx = &tracks[left.track_idx].tx;
    let right_tx = &tracks[right.track_idx].tx;

    left_tx
        .chrom
        .cmp(&right_tx.chrom)
        .then_with(|| strand_rank(left_tx.strand).cmp(&strand_rank(right_tx.strand)))
        .then_with(|| left.interval.start.cmp(&right.interval.start))
        .then_with(|| left.interval.end.cmp(&right.interval.end))
        .then_with(|| left.track_idx.cmp(&right.track_idx))
}

fn same_exon_partition(left: &ExonRecord, right: &ExonRecord, tracks: &[Track]) -> bool {
    let left_tx = &tracks[left.track_idx].tx;
    let right_tx = &tracks[right.track_idx].tx;
    left_tx.chrom == right_tx.chrom && left_tx.strand == right_tx.strand
}

fn exon_overlap_pair_candidates(tracks: &[Track]) -> Vec<PairOverlap> {
    let mut exons: Vec<ExonRecord> = Vec::new();
    for (track_idx, track) in tracks.iter().enumerate() {
        for &interval in &track.tx.exons {
            if !interval.is_empty() {
                exons.push(ExonRecord {
                    track_idx,
                    interval,
                });
            }
        }
    }

    exons.sort_by(|left, right| exon_record_cmp(left, right, tracks));

    let mut active: Vec<usize> = Vec::new();
    let mut overlaps: Vec<PairOverlap> = Vec::new();

    for current_idx in 0..exons.len() {
        let current = exons[current_idx];
        active.retain(|&active_idx| {
            let previous = exons[active_idx];
            same_exon_partition(&previous, &current, tracks)
                && previous.interval.end > current.interval.start
        });

        for &active_idx in &active {
            let previous = exons[active_idx];
            if previous.track_idx == current.track_idx {
                continue;
            }

            let overlap = previous.interval.overlap_len(current.interval);
            if overlap == 0 {
                continue;
            }

            let pair = TrackPair::new(previous.track_idx, current.track_idx);
            overlaps.push(PairOverlap {
                pair,
                exon_overlap: overlap,
            });
        }

        active.push(current_idx);
    }

    overlaps.sort_by(|left, right| left.pair.cmp(&right.pair));

    let mut pairs: Vec<PairOverlap> = Vec::with_capacity(overlaps.len());
    for overlap in overlaps {
        if let Some(last) = pairs.last_mut() {
            if last.pair == overlap.pair {
                last.exon_overlap += overlap.exon_overlap;
                continue;
            }
        }
        pairs.push(overlap);
    }

    pairs
}

fn can_skip_no_exon_overlap(cutoff: f64, intron_weight: f64) -> bool {
    if !cutoff.is_finite() || !intron_weight.is_finite() || intron_weight < 0.0 {
        return false;
    }

    cutoff <= 1.0 / (1.0 + intron_weight)
}

fn should_use_sparse_pair_candidates(tracks_len: usize, cutoff: f64, intron_weight: f64) -> bool {
    tracks_len >= SPARSE_PAIR_THRESHOLD && can_skip_no_exon_overlap(cutoff, intron_weight)
}

fn merge_subreads(src: usize, dst: usize, tracks: &mut [Track]) {
    let subs: Vec<ReadInstance> = tracks[src].subreads.iter().cloned().collect();
    tracks[dst].subreads.extend(subs);
}

fn merge_tracks_by_name(tracks: Vec<Track>) -> Vec<Track> {
    let mut merged: Vec<Track> = Vec::with_capacity(tracks.len());
    let mut positions: HashMap<(TrackSource, String, String), usize> = HashMap::new();

    for track in tracks {
        // A read label identifies a molecule, not a unique alignment. Keep
        // structurally distinct alignments separate during clustering; only
        // exact structural copies of the same source/name are coalesced.
        let key = (
            track.source,
            track.tx.name.clone(),
            crate::identity::novel_isoform_id(&track.tx),
        );
        if let Some(&idx) = positions.get(&key) {
            merged[idx].subreads.extend(track.subreads);
            continue;
        }

        positions.insert(key, merged.len());
        merged.push(track);
    }

    merged
}

fn split_reference_and_read_tracks(tracks: Vec<Track>) -> (Vec<Track>, Vec<Track>) {
    let mut refs: Vec<Track> = Vec::new();
    let mut reads: Vec<Track> = Vec::new();

    for track in tracks {
        if track.is_reference() {
            refs.push(track);
        } else {
            reads.push(track);
        }
    }

    (refs, reads)
}

fn readall(tracks: &[Track]) -> HashSet<usize> {
    let mut instances = HashSet::new();
    for track in tracks {
        for sub in &track.subreads {
            instances.insert(sub.index);
        }
    }
    instances
}

fn readall_subset(tracks: &[Track], keep: &HashSet<usize>) -> HashSet<usize> {
    let mut instances = HashSet::new();
    for &idx in keep {
        let track = &tracks[idx];
        for sub in &track.subreads {
            instances.insert(sub.index);
        }
    }
    instances
}

fn should_drop_read(track: &Track, mode: DistanceMode, sw_score: i64) -> bool {
    mode == DistanceMode::Ratio || sw_score < 0 || i64::from(track.tx.score) < sw_score
}

fn filter_pair(
    tracks: &mut [Track],
    drop: &mut HashSet<usize>,
    pair: PairDistance,
    params: FilterParams,
) {
    let i = pair.left;
    let j = pair.right;
    let distance = combined_distance_from_overlaps(
        &tracks[i],
        &tracks[j],
        params.mode,
        params.intron_weight,
        pair.exon_overlap,
        pair.intron_overlap,
    );
    if distance >= params.cutoff {
        return;
    }

    let li = tracks[i].exon_len;
    let lj = tracks[j].exon_len;

    match (tracks[i].is_reference(), tracks[j].is_reference()) {
        (true, true) => {}
        (true, false) => {
            if should_drop_read(&tracks[j], params.mode, params.sw_score) {
                drop.insert(j);
                merge_subreads(j, i, tracks);
            }
        }
        (false, true) => {
            if should_drop_read(&tracks[i], params.mode, params.sw_score) {
                drop.insert(i);
                merge_subreads(i, j, tracks);
            }
        }
        (false, false) => match li.cmp(&lj) {
            std::cmp::Ordering::Less => {
                if should_drop_read(&tracks[i], params.mode, params.sw_score) {
                    drop.insert(i);
                    merge_subreads(i, j, tracks);
                }
            }
            std::cmp::Ordering::Equal => {
                drop.insert(i);
                merge_subreads(i, j, tracks);
            }
            std::cmp::Ordering::Greater => {
                if should_drop_read(&tracks[j], params.mode, params.sw_score) {
                    drop.insert(j);
                    merge_subreads(j, i, tracks);
                }
            }
        },
    }
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
    let params = FilterParams {
        mode,
        cutoff,
        intron_weight,
        sw_score,
    };

    if should_use_sparse_pair_candidates(tracks.len(), cutoff, intron_weight) {
        for candidate in exon_overlap_pair_candidates(&tracks) {
            let i = candidate.pair.left;
            let j = candidate.pair.right;
            let intron_overlap = overlap_bp(&tracks[i].introns, &tracks[j].introns);
            filter_pair(
                &mut tracks,
                &mut drop,
                PairDistance {
                    left: i,
                    right: j,
                    exon_overlap: candidate.exon_overlap,
                    intron_overlap,
                },
                params,
            );
        }
    } else {
        for i in 0..tracks.len() {
            for j in (i + 1)..tracks.len() {
                let exon_overlap = exonic_overlap_bp(&tracks[i].tx, &tracks[j].tx);
                let intron_overlap = overlap_bp(&tracks[i].introns, &tracks[j].introns);
                filter_pair(
                    &mut tracks,
                    &mut drop,
                    PairDistance {
                        left: i,
                        right: j,
                        exon_overlap,
                        intron_overlap,
                    },
                    params,
                );
            }
        }
    }

    let mut keep: HashSet<usize> = (0..tracks.len()).collect();
    for idx in &drop {
        keep.remove(idx);
    }

    for (idx, track) in tracks.iter().enumerate() {
        if track.is_reference() {
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
        let mut pos: HashMap<usize, usize> = HashMap::new();
        for (idx, track) in tracks.iter().enumerate() {
            if track.is_read() {
                for instance in &track.subreads {
                    pos.insert(instance.index, idx);
                }
            }
        }
        for instance in missed {
            if let Some(idx) = pos.get(&instance).copied() {
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

fn build_read_to_isoform(isoforms: &[Track]) -> Vec<(String, String)> {
    let mut pairs: Vec<(String, String)> = Vec::new();
    for track in isoforms {
        for subread in &track.subreads {
            pairs.push((subread.name.clone(), track.tx.name.clone()));
        }
    }
    pairs.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));
    pairs
}

fn update_name2(isoforms: &mut [Track], mode: Name2Mode) {
    if mode == Name2Mode::None {
        for track in isoforms.iter_mut() {
            track.tx.metadata_mut().set_name2("none");
        }
        return;
    }

    let values: Vec<String> = {
        let mut occurrence: HashMap<usize, u32> = HashMap::new();
        for track in isoforms.iter() {
            for subread in &track.subreads {
                *occurrence.entry(subread.index).or_insert(0) += 1;
            }
        }

        isoforms
            .iter()
            .map(|track| {
                let mut subreads: Vec<&ReadInstance> = track.subreads.iter().collect();
                subreads.sort_unstable_by(|left, right| {
                    left.name
                        .cmp(&right.name)
                        .then_with(|| left.index.cmp(&right.index))
                });
                let mut coverage = 0.0f64;
                for subread in &subreads {
                    let denom = occurrence.get(&subread.index).copied().unwrap_or(0);
                    if denom > 0 {
                        coverage += 1.0f64 / denom as f64;
                    }
                }

                match mode {
                    Name2Mode::Full => crate::identity::encode_name2(
                        subreads.iter().map(|subread| subread.name.as_str()),
                        coverage,
                    )
                    .expect("read IDs were validated before clustering"),
                    Name2Mode::Coverage => format!("|{coverage}"),
                    Name2Mode::None => unreachable!("handled above"),
                }
            })
            .collect()
    };

    for (track, value) in isoforms.iter_mut().zip(values) {
        track.tx.metadata_mut().set_name2(value);
    }
}

struct PartitionResult {
    isoforms: Vec<Transcript>,
    pairs: Vec<(String, String)>,
    represented_read_indices: HashSet<usize>,
    unused_read_indices: Vec<usize>,
}

struct WorkItem {
    index: usize,
    ref_indices: Vec<usize>,
    read_indices: Vec<usize>,
}

fn strand_rank(strand: Strand) -> u8 {
    match strand {
        Strand::Plus => 0,
        Strand::Minus => 1,
        Strand::Unknown => 2,
    }
}

fn sort_tracks_by_coord(tracks: &mut [Track]) {
    tracks.sort_by(|left, right| {
        left.tx
            .chrom
            .cmp(&right.tx.chrom)
            .then_with(|| left.tx.tx_start.cmp(&right.tx.tx_start))
            .then_with(|| left.tx.tx_end.cmp(&right.tx.tx_end))
            .then_with(|| strand_rank(left.tx.strand).cmp(&strand_rank(right.tx.strand)))
            .then_with(|| left.tx.exons.cmp(&right.tx.exons))
            .then_with(|| left.source.cmp(&right.source))
            .then_with(|| left.tx.name.cmp(&right.tx.name))
    });
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
    let read_count = tracks.iter().filter(|track| track.is_read()).count();
    if read_count <= batch_size {
        return cluster_once(tracks, options);
    }

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
    options: ClusterOptions,
) -> PartitionResult {
    let mut records: Vec<Track> = Vec::with_capacity(ref_indices.len() + read_indices.len());
    for &idx in ref_indices {
        records.push(Track::reference(references[idx].clone()));
    }
    for &idx in read_indices {
        records.push(Track::read(reads[idx].clone(), idx));
    }

    sort_tracks_by_coord(&mut records);
    let locus_records: Vec<Transcript> = records.iter().map(|track| track.tx.clone()).collect();
    let loci = cluster_by_span(&locus_records, StrandMode::Match);
    let mut records: Vec<Option<Track>> = records.into_iter().map(Some).collect();

    let mut isoforms: Vec<Transcript> = Vec::new();
    let mut pairs: Vec<(String, String)> = Vec::new();
    let mut represented_read_indices = HashSet::new();
    let mut unused_read_indices = Vec::new();

    for locus in loci {
        let mut tracks: Vec<Track> = Vec::with_capacity(locus.members.len());
        for idx in locus.members {
            tracks.push(records[idx].take().expect("record already consumed"));
        }

        if !tracks.iter().any(Track::is_reference) {
            unused_read_indices.extend(
                tracks
                    .iter()
                    .flat_map(|track| track.subreads.iter().map(|subread| subread.index)),
            );
            continue;
        }

        let mut tracks = batch_overlap_merge(tracks, options);

        for track in &mut tracks {
            if track.is_read() {
                track.tx.name = crate::identity::novel_isoform_id(&track.tx);
            }
        }
        let mut tracks = merge_tracks_by_name(tracks);
        tracks.sort_by(|left, right| left.tx.name.cmp(&right.tx.name));

        update_name2(&mut tracks, options.name2_mode);
        represented_read_indices.extend(
            tracks
                .iter()
                .flat_map(|track| track.subreads.iter().map(|subread| subread.index)),
        );
        pairs.extend(build_read_to_isoform(&tracks));
        isoforms.extend(tracks.into_iter().map(|track| track.tx));
    }

    PartitionResult {
        isoforms,
        pairs,
        represented_read_indices,
        unused_read_indices,
    }
}

pub fn cluster(
    reads: &[Transcript],
    references: Option<&[Transcript]>,
    threads: usize,
) -> ClusterResult {
    try_cluster(reads, references, threads)
        .unwrap_or_else(|error| panic!("invalid overlap-clustering options: {error}"))
}

/// Cluster with default options, returning invalid configuration errors.
pub fn try_cluster(
    reads: &[Transcript],
    references: Option<&[Transcript]>,
    threads: usize,
) -> Result<ClusterResult, crate::config::ParameterError> {
    try_cluster_with_options(reads, references, threads, ClusterOptions::default())
}

pub fn cluster_with_options(
    reads: &[Transcript],
    references: Option<&[Transcript]>,
    threads: usize,
    options: ClusterOptions,
) -> ClusterResult {
    try_cluster_with_options(reads, references, threads, options)
        .unwrap_or_else(|error| panic!("invalid overlap-clustering options: {error}"))
}

/// Cluster with explicit options, returning invalid configuration errors.
pub fn try_cluster_with_options(
    reads: &[Transcript],
    references: Option<&[Transcript]>,
    threads: usize,
    options: ClusterOptions,
) -> Result<ClusterResult, crate::config::ParameterError> {
    let threads = crate::config::WorkerThreads::new(threads)?.get();
    options.validate()?;
    crate::identity::validate_read_ids(reads)
        .map_err(crate::config::ParameterError::invalid_identity)?;
    let references = match references {
        Some(references) => references,
        None => {
            return Ok(ClusterResult {
                isoforms: Vec::new(),
                read_to_isoform: Vec::new(),
                unused: reads.to_vec(),
            });
        }
    };
    crate::identity::validate_reference_ids(references)
        .map_err(crate::config::ParameterError::invalid_identity)?;

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
    let mut unmatched_read_indices = Vec::new();
    for (idx, read) in reads.iter().enumerate() {
        let key = PartitionKey {
            chrom: read.chrom.clone(),
            strand: read.strand,
        };
        if !refs_by_key.contains_key(&key) {
            unmatched_read_indices.push(idx);
            continue;
        }
        reads_by_key.entry(key).or_default().push(idx);
    }

    let mut all_isoforms: Vec<Transcript> = Vec::new();
    let mut all_pairs: Vec<(String, String)> = Vec::new();
    let mut represented_read_indices = HashSet::new();
    let mut unused_read_indices = unmatched_read_indices;

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
        represented_read_indices.extend(part.represented_read_indices);
        unused_read_indices.extend(part.unused_read_indices);
    }

    unused_read_indices.sort_unstable();
    let unused_instance_set: HashSet<usize> = unused_read_indices.iter().copied().collect();
    assert_eq!(
        unused_instance_set.len(),
        unused_read_indices.len(),
        "overlap clustering classified a read instance as unused more than once"
    );
    assert!(
        represented_read_indices.is_disjoint(&unused_instance_set),
        "overlap clustering classified a read instance as both represented and unused"
    );
    assert!(
        represented_read_indices
            .iter()
            .chain(unused_instance_set.iter())
            .all(|index| *index < reads.len()),
        "overlap clustering produced an out-of-range read instance"
    );
    assert_eq!(
        represented_read_indices.len() + unused_instance_set.len(),
        reads.len(),
        "overlap clustering violated read-instance conservation"
    );

    all_pairs.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));
    all_isoforms.sort_by(|left, right| left.name.cmp(&right.name));
    crate::identity::validate_isoform_ids(&all_isoforms)
        .map_err(crate::config::ParameterError::invalid_identity)?;
    let mut all_unused: Vec<Transcript> = unused_read_indices
        .into_iter()
        .map(|index| reads[index].clone())
        .collect();
    all_unused.sort_by(crate::identity::transcript_order);

    Ok(ClusterResult {
        isoforms: all_isoforms,
        read_to_isoform: all_pairs,
        unused: all_unused,
    })
}

#[cfg(test)]
mod tests {
    use crate::model::{Bed12Attrs, Coord, Interval, Strand, Transcript};

    use super::*;

    #[test]
    fn name2_coverage_sum_is_independent_of_hash_iteration_order() {
        let expected = (1..=12)
            .map(|denominator| 1.0 / f64::from(denominator))
            .sum::<f64>();
        let expected_payload = format!("|{expected}");

        for _ in 0..64 {
            let mut tracks = (0..12)
                .map(|index| {
                    Track::reference(make_tx(
                        &format!("isoform_{index:02}"),
                        Strand::Plus,
                        &[(100, 200)],
                        "isoform_anno",
                        100,
                    ))
                })
                .collect::<Vec<_>>();
            for read_index in 0..12 {
                let read = ReadInstance {
                    index: read_index,
                    name: format!("read_{read_index:02}"),
                };
                for track in tracks.iter_mut().take(read_index + 1) {
                    track.subreads.insert(read.clone());
                }
            }

            update_name2(&mut tracks, Name2Mode::Coverage);
            assert_eq!(
                tracks[0].tx.metadata().name2_field(),
                Some(expected_payload.as_str())
            );
        }
    }

    fn make_tx(
        name: &str,
        strand: Strand,
        exons: &[(u32, u32)],
        ttype: &str,
        score: u32,
    ) -> Transcript {
        make_tx_on("chr1", name, strand, exons, ttype, score)
    }

    fn make_tx_on(
        chrom: &str,
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
            chrom.to_owned(),
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

    fn make_plain_bed12_tx(
        chrom: &str,
        name: &str,
        strand: Strand,
        exons: &[(u32, u32)],
        score: u32,
    ) -> Transcript {
        let tx_start = exons.iter().map(|(s, _)| *s).min().unwrap_or(0);
        let tx_end = exons.iter().map(|(_, e)| *e).max().unwrap_or(0);
        let exons = exons
            .iter()
            .map(|(s, e)| Interval::new(Coord::new(*s), Coord::new(*e)).unwrap())
            .collect::<Vec<_>>();

        Transcript::new(
            chrom.to_owned(),
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
                extra_fields: Vec::new(),
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
    fn plain_bed12_reference_is_protected_by_source_not_ttype() {
        let refs = vec![make_plain_bed12_tx(
            "chr1",
            "ref",
            Strand::Plus,
            &[(100, 110), (120, 130), (140, 150)],
            100,
        )];
        let reads = vec![make_tx(
            "read",
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

        assert!(result.isoforms.iter().any(|tx| tx.name == "ref"));
        assert!(!result.isoforms.iter().any(|tx| tx.name == "read"));
        assert_eq!(
            result.read_to_isoform,
            vec![("read".to_owned(), "ref".to_owned())]
        );
        assert!(result.unused.is_empty());
    }

    #[test]
    fn unmatched_reads_are_reported_as_unused() {
        let refs = vec![make_plain_bed12_tx(
            "chr1",
            "ref",
            Strand::Plus,
            &[(100, 110), (120, 130)],
            100,
        )];
        let reads = vec![
            make_tx(
                "read_match",
                Strand::Plus,
                &[(100, 110), (120, 130)],
                "nanopore_read",
                0,
            ),
            make_tx(
                "read_disjoint",
                Strand::Plus,
                &[(500, 510), (520, 530)],
                "nanopore_read",
                0,
            ),
            make_tx(
                "read_wrong_strand",
                Strand::Minus,
                &[(100, 110), (120, 130)],
                "nanopore_read",
                0,
            ),
            make_tx_on(
                "chr2",
                "read_wrong_chrom",
                Strand::Plus,
                &[(100, 110), (120, 130)],
                "nanopore_read",
                0,
            ),
        ];

        let result = cluster_with_options(
            &reads,
            Some(&refs),
            1,
            ClusterOptions {
                name2_mode: Name2Mode::Coverage,
                ..ClusterOptions::default()
            },
        );

        assert_eq!(
            result.read_to_isoform,
            vec![("read_match".to_owned(), "ref".to_owned())]
        );
        let unused: HashSet<&str> = result.unused.iter().map(|tx| tx.name.as_str()).collect();
        assert_eq!(unused.len(), 3);
        assert!(unused.contains("read_disjoint"));
        assert!(unused.contains("read_wrong_strand"));
        assert!(unused.contains("read_wrong_chrom"));
    }

    #[test]
    fn duplicate_read_ids_are_conserved_as_distinct_instances() {
        let refs = vec![make_plain_bed12_tx(
            "chr1",
            "ref",
            Strand::Plus,
            &[(100, 110), (120, 130)],
            100,
        )];
        let reads = vec![
            make_plain_bed12_tx(
                "chr1",
                "duplicate_id",
                Strand::Plus,
                &[(100, 110), (120, 130)],
                1,
            ),
            make_plain_bed12_tx(
                "chr2",
                "duplicate_id",
                Strand::Plus,
                &[(100, 110), (120, 130)],
                1,
            ),
        ];

        let result = cluster_with_options(
            &reads,
            Some(&refs),
            1,
            ClusterOptions {
                name2_mode: Name2Mode::Coverage,
                ..ClusterOptions::default()
            },
        );

        assert_eq!(
            result.read_to_isoform,
            vec![("duplicate_id".to_owned(), "ref".to_owned())]
        );
        assert_eq!(result.unused.len(), 1);
        assert_eq!(result.unused[0].chrom, "chr2");
        assert_eq!(result.unused[0].name, "duplicate_id");
    }

    #[test]
    fn sw_score_minus_one_keeps_ratio_short_truncation_merge() {
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
        let no_signal_merge = cluster_with_options(
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
        let no_signal_merge_reads: HashSet<&str> = no_signal_merge
            .read_to_isoform
            .iter()
            .map(|(read_id, _)| read_id.as_str())
            .collect();
        assert_eq!(merged_reads.len(), 2);
        assert_eq!(no_signal_merge_reads.len(), 2);

        let no_signal_merge_targets: HashSet<&str> = no_signal_merge
            .read_to_isoform
            .iter()
            .map(|(_, isoform_id)| isoform_id.as_str())
            .collect();

        assert_eq!(no_signal_merge_targets.len(), 1);

        let merged_short_targets: HashSet<&str> = merged
            .read_to_isoform
            .iter()
            .filter(|(read_id, _)| read_id == "read_short")
            .map(|(_, isoform_id)| isoform_id.as_str())
            .collect();
        let no_signal_merge_short_targets: HashSet<&str> = no_signal_merge
            .read_to_isoform
            .iter()
            .filter(|(read_id, _)| read_id == "read_short")
            .map(|(_, isoform_id)| isoform_id.as_str())
            .collect();

        assert_eq!(merged_short_targets.len(), 1);
        assert!(merged_short_targets.contains("ref"));
        assert_eq!(no_signal_merge_short_targets.len(), 1);
        assert!(no_signal_merge_short_targets.contains("ref"));
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
        assert!(sl_targets
            .iter()
            .all(|id| id.starts_with(crate::identity::NOVEL_ISOFORM_PREFIX)));
    }

    #[test]
    fn sw_score_above_cutoff_keeps_sl_read_as_its_own_track() {
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
                DEFAULT_SW_SCORE as u32 + 1,
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
        assert!(sl_targets
            .iter()
            .all(|id| id.starts_with(crate::identity::NOVEL_ISOFORM_PREFIX)));
    }

    #[test]
    fn sw_score_minus_one_removes_overlap_sl_protection() {
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
                DEFAULT_SW_SCORE as u32 + 1,
            ),
        ];

        let protected = cluster_with_options(&reads, Some(&refs), 1, ClusterOptions::default());
        let protected_targets: HashSet<&str> = protected
            .read_to_isoform
            .iter()
            .filter(|(read_id, _)| read_id == "read_sl")
            .map(|(_, isoform_id)| isoform_id.as_str())
            .collect();
        assert_eq!(protected_targets.len(), 1);
        assert!(protected_targets
            .iter()
            .all(|id| id.starts_with(crate::identity::NOVEL_ISOFORM_PREFIX)));

        for batch_size in [0, 1] {
            let no_signal = cluster_with_options(
                &reads,
                Some(&refs),
                1,
                ClusterOptions {
                    sw_score: -1,
                    batch_size,
                    ..ClusterOptions::default()
                },
            );
            let no_signal_targets: HashSet<&str> = no_signal
                .read_to_isoform
                .iter()
                .filter(|(read_id, _)| read_id == "read_sl")
                .map(|(_, isoform_id)| isoform_id.as_str())
                .collect();

            assert_eq!(no_signal_targets.len(), 1);
            assert!(no_signal_targets.contains("ref"));
        }
    }

    #[test]
    fn exon_candidate_generation_is_sparse_and_ordered() {
        let tracks = vec![
            Track::reference(make_tx(
                "ref",
                Strand::Plus,
                &[(100, 110), (120, 130)],
                "isoform_anno",
                100,
            )),
            Track::read(
                make_tx(
                    "read_overlap",
                    Strand::Plus,
                    &[(105, 115), (125, 135)],
                    "nanopore_read",
                    0,
                ),
                0,
            ),
            Track::read(
                make_tx(
                    "read_disjoint",
                    Strand::Plus,
                    &[(200, 210), (220, 230)],
                    "nanopore_read",
                    0,
                ),
                1,
            ),
        ];

        let candidates = exon_overlap_pair_candidates(&tracks);

        assert_eq!(candidates.len(), 1);
        assert_eq!(candidates[0].pair, TrackPair { left: 0, right: 1 });
        assert_eq!(candidates[0].exon_overlap, 10);
    }

    #[test]
    fn high_cutoff_preserves_no_exon_overlap_intron_merge_behavior() {
        let refs = vec![make_tx(
            "ref",
            Strand::Plus,
            &[(100, 110), (130, 140)],
            "isoform_anno",
            100,
        )];
        let reads = vec![make_tx(
            "read_no_exon_overlap",
            Strand::Plus,
            &[(111, 120), (140, 150)],
            "nanopore_read",
            0,
        )];

        let result = cluster_with_options(
            &reads,
            Some(&refs),
            1,
            ClusterOptions {
                cutoff1: 0.95,
                name2_mode: Name2Mode::Coverage,
                ..ClusterOptions::default()
            },
        );

        assert_eq!(
            result.read_to_isoform,
            vec![("read_no_exon_overlap".to_owned(), "ref".to_owned())]
        );
    }

    #[test]
    fn overlap_cluster_is_deterministic_across_thread_counts() {
        let refs = vec![
            make_tx_on(
                "chr1",
                "ref1",
                Strand::Plus,
                &[(100, 110), (120, 130)],
                "isoform_anno",
                100,
            ),
            make_tx_on(
                "chr2",
                "ref2",
                Strand::Plus,
                &[(100, 110), (120, 130)],
                "isoform_anno",
                100,
            ),
        ];
        let reads = vec![
            make_tx_on(
                "chr1",
                "read1",
                Strand::Plus,
                &[(100, 110), (120, 130)],
                "nanopore_read",
                0,
            ),
            make_tx_on(
                "chr2",
                "read2",
                Strand::Plus,
                &[(100, 110), (120, 130)],
                "nanopore_read",
                0,
            ),
        ];

        let single = cluster_with_options(&reads, Some(&refs), 1, ClusterOptions::default());
        let multi = cluster_with_options(&reads, Some(&refs), 4, ClusterOptions::default());

        assert_eq!(single.isoforms, multi.isoforms);
        assert_eq!(single.read_to_isoform, multi.read_to_isoform);
        assert_eq!(single.unused, multi.unused);
    }

    #[test]
    fn stable_novel_ids_and_safe_name2_are_input_order_independent() {
        let refs = vec![make_tx(
            "ref",
            Strand::Plus,
            &[(100, 110), (120, 130), (140, 150)],
            "isoform_anno",
            100,
        )];
        let reads = vec![
            make_tx(
                "read,z|%",
                Strand::Plus,
                &[(120, 130), (140, 150)],
                "nanopore_read",
                DEFAULT_SW_SCORE as u32 + 1,
            ),
            make_tx(
                "read,a",
                Strand::Plus,
                &[(120, 130), (140, 150)],
                "nanopore_read",
                DEFAULT_SW_SCORE as u32 + 1,
            ),
        ];

        let forward = cluster_with_options(&reads, Some(&refs), 1, ClusterOptions::default());
        let mut reversed_reads = reads.clone();
        reversed_reads.reverse();
        let reversed =
            cluster_with_options(&reversed_reads, Some(&refs), 4, ClusterOptions::default());

        assert_eq!(forward.isoforms, reversed.isoforms);
        assert_eq!(forward.read_to_isoform, reversed.read_to_isoform);
        let targets: HashSet<&str> = forward
            .read_to_isoform
            .iter()
            .filter(|(read, _)| read == "read,z|%" || read == "read,a")
            .map(|(_, isoform)| isoform.as_str())
            .collect();
        assert_eq!(targets.len(), 1);
        let novel_id = *targets.iter().next().unwrap();
        assert!(novel_id.starts_with(crate::identity::NOVEL_ISOFORM_PREFIX));
        let novel = forward
            .isoforms
            .iter()
            .find(|tx| tx.name == novel_id)
            .unwrap();
        assert_eq!(
            crate::identity::decode_name2(&novel.extra_fields[0])
                .unwrap()
                .into_iter()
                .collect::<HashSet<_>>(),
            HashSet::from(["read,z|%".to_owned(), "read,a".to_owned()])
        );
    }

    #[test]
    fn overlap_clustering_rejects_duplicate_reference_ids() {
        let refs = vec![
            make_tx("same", Strand::Plus, &[(0, 10)], "isoform_anno", 100),
            make_tx("same", Strand::Plus, &[(20, 30)], "isoform_anno", 100),
        ];
        let error = try_cluster(&[], Some(&refs), 1).unwrap_err();
        assert!(error.to_string().contains("duplicate reference isoform id"));
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
        let oversized_batch = cluster_with_options(
            &reads,
            Some(&refs),
            1,
            ClusterOptions {
                batch_size: 100_000,
                batch_rounds: 10,
                ..ClusterOptions::default()
            },
        );

        assert_eq!(single_pass.read_to_isoform, batched.read_to_isoform);
        assert_eq!(single_pass.read_to_isoform, oversized_batch.read_to_isoform);

        let single_names: Vec<&str> = single_pass
            .isoforms
            .iter()
            .map(|tx| tx.name.as_str())
            .collect();
        let batch_names: Vec<&str> = batched.isoforms.iter().map(|tx| tx.name.as_str()).collect();
        let oversized_names: Vec<&str> = oversized_batch
            .isoforms
            .iter()
            .map(|tx| tx.name.as_str())
            .collect();
        assert_eq!(single_names, batch_names);
        assert_eq!(single_names, oversized_names);
    }

    #[test]
    fn rejects_invalid_options_at_library_boundary() {
        assert!(try_cluster_with_options(&[], None, 0, ClusterOptions::default()).is_err());

        for options in [
            ClusterOptions {
                cutoff1: f64::NAN,
                ..ClusterOptions::default()
            },
            ClusterOptions {
                cutoff2: 1.1,
                ..ClusterOptions::default()
            },
            ClusterOptions {
                intron_weight: -0.1,
                ..ClusterOptions::default()
            },
            ClusterOptions {
                intron_weight: f64::INFINITY,
                ..ClusterOptions::default()
            },
        ] {
            assert!(try_cluster_with_options(&[], None, 1, options).is_err());
        }
    }
}
