use std::collections::{HashMap, HashSet};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use crate::cluster::result::ClusterResult;
use crate::interval::{cluster_by_span, StrandMode};
use crate::model::{Coord, Interval, Strand, Transcript};

#[derive(Clone, Debug)]
struct Track {
    tx: Transcript,
    source: TrackSource,
    subreads: HashSet<ReadInstance>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
enum TrackSource {
    Reference,
    Read,
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct ReadInstance {
    index: usize,
    name: String,
}

impl Track {
    fn reference(tx: Transcript) -> Self {
        Self {
            tx,
            source: TrackSource::Reference,
            subreads: HashSet::new(),
        }
    }

    fn read(tx: Transcript, index: usize) -> Self {
        let subreads = HashSet::from([ReadInstance {
            index,
            name: tx.name.clone(),
        }]);
        Self {
            tx,
            source: TrackSource::Read,
            subreads,
        }
    }

    fn is_reference(&self) -> bool {
        self.source == TrackSource::Reference
    }

    fn is_read(&self) -> bool {
        self.source == TrackSource::Read
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct PartitionKey {
    chrom: String,
    strand: Strand,
}

pub const DEFAULT_SW_SCORE: i64 = -1;
pub const DEFAULT_JUNCTION_CORRECTION_MIN_SUPPORT: u32 = 5;
pub const DEFAULT_JUNCTION_CORRECTION_OFFSET: u32 = 10;
pub const DEFAULT_SL_PARTIAL_FIVE_PRIME_END_OFFSET: u32 = 15;
pub const DEFAULT_SL_SAME_JUNCTION_FIVE_PRIME_END_OFFSET: u32 = 25;
pub const DEFAULT_SL_FIVE_PRIME_CLUSTER_OFFSET: u32 = 15;
pub const DEFAULT_MIN_SL_FIVE_PRIME_CLUSTER_SUPPORT: usize = 2;
pub const DEFAULT_SAME_JUNCTION_THREE_PRIME_END_OFFSET: u32 = 50;
pub const DEFAULT_MIN_THREE_PRIME_CLUSTER_SUPPORT: usize = 5;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct JunctionClusterSummary {
    pub input_reads: usize,
    pub represented_reads: usize,
    pub mapping_rows: usize,
    pub rare_reads: usize,
    pub unmatched_reads: usize,
    pub unused_reads: usize,
}

impl JunctionClusterSummary {
    pub(crate) fn emit(self) {
        eprintln!(
            "clusterj: input_reads={} represented_reads={} mapping_rows={} rare_reads={} unmatched_reads={} unused_reads={}",
            self.input_reads,
            self.represented_reads,
            self.mapping_rows,
            self.rare_reads,
            self.unmatched_reads,
            self.unused_reads
        );
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct JunctionCorrectionOptions {
    pub min_support: u32,
    pub offset: u32,
}

impl Default for JunctionCorrectionOptions {
    fn default() -> Self {
        Self {
            min_support: DEFAULT_JUNCTION_CORRECTION_MIN_SUPPORT,
            offset: DEFAULT_JUNCTION_CORRECTION_OFFSET,
        }
    }
}

impl JunctionCorrectionOptions {
    /// Validate junction-correction support and offset domains.
    pub fn validate(self) -> Result<(), crate::config::ParameterError> {
        crate::config::WeightedMinimumSupport::new(
            "junction correction minimum support",
            self.min_support,
        )?;
        let _ = crate::config::BasePairOffset::new(self.offset);
        Ok(())
    }
}

/// Runtime bounds for one `clusterj` invocation (downsample, heartbeat).
///
/// Library callers default to no per-locus cap and no heartbeat. The
/// standalone CLI enables the same 5,000-read locus cap used by `flow`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ClusterjRuntimeOptions {
    /// Reservoir cap per overlapping locus; zero disables downsampling.
    pub max_reads_per_locus: usize,
    /// Base seed mixed with chrom, strand, and locus span.
    pub downsample_seed: u64,
    /// Heartbeat interval in seconds; zero disables.
    pub heartbeat_seconds: u64,
    /// How many in-flight partitions to print when a heartbeat sees no progress.
    pub heartbeat_top: usize,
}

impl Default for ClusterjRuntimeOptions {
    fn default() -> Self {
        Self {
            max_reads_per_locus: 0,
            downsample_seed: 1,
            heartbeat_seconds: 0,
            heartbeat_top: 5,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SlMergeOptions {
    pub partial_five_prime_end_offset: u32,
    pub same_junction_five_prime_end_offset: u32,
    pub five_prime_cluster_offset: u32,
    pub min_five_prime_cluster_support: usize,
}

impl Default for SlMergeOptions {
    fn default() -> Self {
        Self {
            partial_five_prime_end_offset: DEFAULT_SL_PARTIAL_FIVE_PRIME_END_OFFSET,
            same_junction_five_prime_end_offset: DEFAULT_SL_SAME_JUNCTION_FIVE_PRIME_END_OFFSET,
            five_prime_cluster_offset: DEFAULT_SL_FIVE_PRIME_CLUSTER_OFFSET,
            min_five_prime_cluster_support: DEFAULT_MIN_SL_FIVE_PRIME_CLUSTER_SUPPORT,
        }
    }
}

impl SlMergeOptions {
    /// Validate SL offset and support domains.
    pub fn validate(self) -> Result<(), crate::config::ParameterError> {
        let _ = crate::config::BasePairOffset::new(self.partial_five_prime_end_offset);
        let _ = crate::config::BasePairOffset::new(self.same_junction_five_prime_end_offset);
        let _ = crate::config::BasePairOffset::new(self.five_prime_cluster_offset);
        crate::config::MinimumSupport::new(
            "SL 5-prime minimum support",
            self.min_five_prime_cluster_support,
        )?;
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ThreePrimeMergeOptions {
    pub same_junction_three_prime_end_offset: u32,
    pub three_prime_cluster_offset: u32,
    pub min_three_prime_cluster_support: usize,
}

impl ThreePrimeMergeOptions {
    fn with_junction_offset(junction_offset: u32) -> Self {
        Self {
            same_junction_three_prime_end_offset: DEFAULT_SAME_JUNCTION_THREE_PRIME_END_OFFSET,
            three_prime_cluster_offset: junction_offset,
            min_three_prime_cluster_support: DEFAULT_MIN_THREE_PRIME_CLUSTER_SUPPORT,
        }
    }

    /// Validate 3-prime offset and support domains.
    pub fn validate(self) -> Result<(), crate::config::ParameterError> {
        let _ = crate::config::BasePairOffset::new(self.same_junction_three_prime_end_offset);
        let _ = crate::config::BasePairOffset::new(self.three_prime_cluster_offset);
        crate::config::MinimumSupport::new(
            "3-prime minimum support",
            self.min_three_prime_cluster_support,
        )?;
        Ok(())
    }
}

impl Default for ThreePrimeMergeOptions {
    fn default() -> Self {
        Self::with_junction_offset(DEFAULT_JUNCTION_CORRECTION_OFFSET)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ResolvedPlatformOptions {
    pub junction_correction: JunctionCorrectionOptions,
    pub sl_options: SlMergeOptions,
    pub three_prime_options: ThreePrimeMergeOptions,
}

impl ResolvedPlatformOptions {
    /// Validate every resolved scientific option.
    pub fn validate(self) -> Result<(), crate::config::ParameterError> {
        self.junction_correction.validate()?;
        self.sl_options.validate()?;
        self.three_prime_options.validate()?;
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum PlatformPreset {
    #[default]
    Generic,
    Rna002,
    Rna004,
}

impl PlatformPreset {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Generic => "generic",
            Self::Rna002 => "rna002",
            Self::Rna004 => "rna004",
        }
    }

    pub fn options(self) -> ResolvedPlatformOptions {
        match self {
            Self::Generic | Self::Rna004 => {
                let junction_correction = JunctionCorrectionOptions::default();
                ResolvedPlatformOptions {
                    junction_correction,
                    sl_options: SlMergeOptions::default(),
                    three_prime_options: ThreePrimeMergeOptions::with_junction_offset(
                        junction_correction.offset,
                    ),
                }
            }
            Self::Rna002 => {
                let junction_correction = JunctionCorrectionOptions {
                    min_support: DEFAULT_JUNCTION_CORRECTION_MIN_SUPPORT,
                    offset: 15,
                };
                ResolvedPlatformOptions {
                    junction_correction,
                    sl_options: SlMergeOptions {
                        partial_five_prime_end_offset: 20,
                        same_junction_five_prime_end_offset:
                            DEFAULT_SL_SAME_JUNCTION_FIVE_PRIME_END_OFFSET,
                        five_prime_cluster_offset: 20,
                        min_five_prime_cluster_support: DEFAULT_MIN_SL_FIVE_PRIME_CLUSTER_SUPPORT,
                    },
                    three_prime_options: ThreePrimeMergeOptions::with_junction_offset(
                        junction_correction.offset,
                    ),
                }
            }
        }
    }
}

impl std::fmt::Display for PlatformPreset {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

impl std::str::FromStr for PlatformPreset {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.eq_ignore_ascii_case("generic") {
            return Ok(Self::Generic);
        }
        if s.eq_ignore_ascii_case("rna002") {
            return Ok(Self::Rna002);
        }
        if s.eq_ignore_ascii_case("rna004") {
            return Ok(Self::Rna004);
        }
        Err(format!(
            "invalid platform preset {s:?}; expected one of: generic, rna002, rna004"
        ))
    }
}

#[allow(clippy::too_many_arguments)]
pub fn resolve_platform_options(
    platform_preset: PlatformPreset,
    junction_correction_offset: Option<u32>,
    junction_correction_min_support: Option<u32>,
    sl_partial_five_prime_offset: Option<u32>,
    sl_same_junction_five_prime_offset: Option<u32>,
    sl_five_prime_cluster_offset: Option<u32>,
    sl_five_prime_min_support: Option<usize>,
    same_junction_three_prime_offset: Option<u32>,
    three_prime_cluster_offset: Option<u32>,
    three_prime_min_support: Option<usize>,
) -> ResolvedPlatformOptions {
    let mut options = platform_preset.options();
    if let Some(offset) = junction_correction_offset {
        options.junction_correction.offset = offset;
    }
    if let Some(min_support) = junction_correction_min_support {
        options.junction_correction.min_support = min_support;
    }
    if let Some(offset) = sl_partial_five_prime_offset {
        options.sl_options.partial_five_prime_end_offset = offset;
    }
    if let Some(offset) = sl_same_junction_five_prime_offset {
        options.sl_options.same_junction_five_prime_end_offset = offset;
    }
    if let Some(offset) = sl_five_prime_cluster_offset {
        options.sl_options.five_prime_cluster_offset = offset;
    }
    if let Some(min_support) = sl_five_prime_min_support {
        options.sl_options.min_five_prime_cluster_support = min_support;
    }
    if three_prime_cluster_offset.is_none() {
        options.three_prime_options.three_prime_cluster_offset = options.junction_correction.offset;
    }
    if let Some(offset) = same_junction_three_prime_offset {
        options
            .three_prime_options
            .same_junction_three_prime_end_offset = offset;
    }
    if let Some(offset) = three_prime_cluster_offset {
        options.three_prime_options.three_prime_cluster_offset = offset;
    }
    if let Some(min_support) = three_prime_min_support {
        options.three_prime_options.min_three_prime_cluster_support = min_support;
    }
    options
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum MergeKind {
    SameJunction,
    FivePrimeTruncation,
    SingleExonContained,
    SingleExonSameFivePrime,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum Name2Mode {
    #[default]
    Full,
    Coverage,
    None,
}

impl Name2Mode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Full => "full",
            Self::Coverage => "coverage",
            Self::None => "none",
        }
    }
}

impl std::fmt::Display for Name2Mode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

impl std::str::FromStr for Name2Mode {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.eq_ignore_ascii_case("full") {
            return Ok(Self::Full);
        }
        if s.eq_ignore_ascii_case("coverage") {
            return Ok(Self::Coverage);
        }
        if s.eq_ignore_ascii_case("none") {
            return Ok(Self::None);
        }
        Err(format!(
            "invalid name2 mode {s:?}; expected one of: full, coverage, none"
        ))
    }
}

fn track_weight(track: &Track, ref_weight: u32, read_weight: u32) -> u32 {
    match track.source {
        TrackSource::Reference => ref_weight,
        TrackSource::Read => read_weight,
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
) -> Option<Vec<Interval>> {
    let start = tx_start.get();
    let end = tx_end.get();

    let mut junctions = junctions.to_vec();
    junctions.sort_unstable();

    if start >= end
        || junctions.is_empty()
        || !junctions.len().is_multiple_of(2)
        || junctions.windows(2).any(|pair| pair[0] >= pair[1])
        || junctions
            .iter()
            .any(|boundary| *boundary <= start || *boundary >= end)
    {
        return None;
    }

    let mut exons: Vec<Interval> = Vec::new();
    exons.push(Interval::new(Coord::new(start), Coord::new(junctions[0])).ok()?);

    let mut idx = 1usize;
    while idx + 1 < junctions.len() {
        let exon =
            Interval::new(Coord::new(junctions[idx]), Coord::new(junctions[idx + 1])).ok()?;
        if exon.is_empty() {
            return None;
        }
        exons.push(exon);
        idx += 2;
    }

    exons.push(Interval::new(Coord::new(*junctions.last()?), Coord::new(end)).ok()?);

    exons.iter().all(|exon| !exon.is_empty()).then_some(exons)
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
        let weight = track_weight(track, 5, 1);
        for &pos in junctions {
            *site_cov.entry(pos).or_insert(0) += weight;
        }
    }

    let (w_to_r, w_to_no) = build_corrected_site_maps(&site_cov, cov_cutoff, pos_cutoff);

    let mut corrected: Vec<Track> = Vec::new();
    let mut rare: Vec<Track> = Vec::new();

    for (idx, mut track) in tracks.into_iter().enumerate() {
        let junctions = &junctions_cache[idx];
        // References are anchors supplied by the reference input. Their metadata is biological
        // annotation only, so never correct or discard them based on transcript-type fields.
        if track.is_reference() {
            corrected.push(track);
            continue;
        }

        if junctions.iter().any(|pos| w_to_no.contains(pos)) {
            rare.push(track);
            continue;
        }

        if !junctions.is_empty() {
            let mut corrected_junctions: Vec<u32> = Vec::with_capacity(junctions.len());
            let mut changed = false;
            for &pos in junctions {
                let corrected_pos = w_to_r.get(&pos).copied().unwrap_or(pos);
                if corrected_pos != pos {
                    changed = true;
                }
                corrected_junctions.push(corrected_pos);
            }

            if changed {
                let Some(corrected_exons) = rebuild_exons_from_junctions(
                    track.tx.tx_start,
                    track.tx.tx_end,
                    &corrected_junctions,
                ) else {
                    rare.push(track);
                    continue;
                };
                if track.tx.metadata().transcript_type().is_some() {
                    track
                        .tx
                        .metadata_mut()
                        .set_transcript_type("nanopore_read_corrected");
                }
                track.tx.exons = corrected_exons;
            }
        }

        corrected.push(track);
    }

    (corrected, rare)
}

fn ordered_boundary_matches(a: &[u32], b: &[u32], offset: u32) -> Vec<(usize, usize)> {
    crate::matching::ordered_one_to_one_matches_by(a.len(), b.len(), |a_idx, b_idx| {
        let delta = a[a_idx].abs_diff(b[b_idx]);
        (delta <= offset).then_some(u64::from(delta))
    })
}

fn junctions_equal(a: &[u32], b: &[u32], offset: u32) -> bool {
    a.len() == b.len()
        && a.iter()
            .zip(b.iter())
            .all(|(left, right)| left.abs_diff(*right) <= offset)
}

fn compare_ei_by_boundary(a: &[u32], reference: &[u32], offset: u32) -> (Vec<usize>, Vec<usize>) {
    let matches = ordered_boundary_matches(a, reference, offset);
    let mut matched_a = vec![false; a.len()];
    let mut matched_reference = vec![false; reference.len()];
    for (a_idx, reference_idx) in matches {
        matched_a[a_idx] = true;
        matched_reference[reference_idx] = true;
    }

    let missed_order = matched_reference
        .iter()
        .enumerate()
        .filter_map(|(idx, matched)| (!matched).then_some(idx))
        .collect();
    let extra_order = matched_a
        .iter()
        .enumerate()
        .filter_map(|(idx, matched)| (!matched).then_some(idx))
        .collect();

    (missed_order, extra_order)
}

fn is_junction_5primer(missed_order: &[usize]) -> bool {
    if missed_order.is_empty() || missed_order[0] != 0 {
        return false;
    }

    let groups = group_consecutive_indices(missed_order);
    groups.len() == 1
}

fn junction_merge_kind(
    short_junctions: &[u32],
    long_junctions: &[u32],
    same_junction_offset: u32,
) -> Option<MergeKind> {
    if short_junctions.is_empty() {
        return None;
    }

    if junctions_equal(short_junctions, long_junctions, same_junction_offset) {
        return Some(MergeKind::SameJunction);
    }

    let (missed_order, extra_order) = compare_ei_by_boundary(short_junctions, long_junctions, 0);
    if missed_order.is_empty() {
        return None;
    }
    (is_junction_5primer(&missed_order) && extra_order.is_empty())
        .then_some(MergeKind::FivePrimeTruncation)
}

fn is_single_exon_in(single: &Transcript, other: &Transcript, other_junctions: &[u32]) -> bool {
    if other_junctions.is_empty() {
        if single.tx_start == other.tx_start && single.tx_end == other.tx_end {
            return true;
        }
        match single.strand {
            Strand::Plus | Strand::Minus => {
                if single.tx_start.get() <= other.tx_start.get()
                    || single.tx_end.get() >= other.tx_end.get()
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
                if single.tx_start.get() < last_junction || single.tx_end.get() > other.tx_end.get()
                {
                    return false;
                }
                true
            }
            Strand::Minus => {
                if single.tx_start.get() < other.tx_start.get()
                    || single.tx_end.get() > last_junction
                {
                    return false;
                }
                true
            }
            Strand::Unknown => false,
        }
    }
}

fn five_prime_end_delta(a: &Transcript, b: &Transcript) -> Option<u32> {
    if a.chrom != b.chrom || a.strand != b.strand {
        return None;
    }

    match a.strand {
        Strand::Plus => Some(a.tx_start.get().abs_diff(b.tx_start.get())),
        Strand::Minus => Some(a.tx_end.get().abs_diff(b.tx_end.get())),
        Strand::Unknown => None,
    }
}

fn five_prime_position(tx: &Transcript) -> Option<u32> {
    match tx.strand {
        Strand::Plus => Some(tx.tx_start.get()),
        Strand::Minus => Some(tx.tx_end.get()),
        Strand::Unknown => None,
    }
}

fn five_prime_ends_match(a: &Transcript, b: &Transcript, offset: u32) -> bool {
    five_prime_end_delta(a, b).is_some_and(|delta| delta <= offset)
}

fn three_prime_end_delta(a: &Transcript, b: &Transcript) -> Option<u32> {
    if a.chrom != b.chrom || a.strand != b.strand {
        return None;
    }

    match a.strand {
        Strand::Plus => Some(a.tx_end.get().abs_diff(b.tx_end.get())),
        Strand::Minus => Some(a.tx_start.get().abs_diff(b.tx_start.get())),
        Strand::Unknown => None,
    }
}

fn three_prime_position(tx: &Transcript) -> Option<u32> {
    match tx.strand {
        Strand::Plus => Some(tx.tx_end.get()),
        Strand::Minus => Some(tx.tx_start.get()),
        Strand::Unknown => None,
    }
}

fn is_sl_supported_read(track: &Track, sw_score: i64) -> bool {
    sw_score >= 0 && track.is_read() && i64::from(track.tx.score) > sw_score
}

fn track_read_support(track: &Track) -> usize {
    if track.is_reference() {
        0
    } else {
        track.subreads.len().max(1)
    }
}

fn build_sl_five_prime_cluster_support(
    tracks: &[Track],
    junctions_cache: &[Vec<u32>],
    sw_score: i64,
    sl_options: SlMergeOptions,
) -> Vec<usize> {
    let mut support = vec![0usize; tracks.len()];
    if sw_score < 0 {
        return support;
    }

    #[derive(Clone, Copy, Debug)]
    struct SlSupportEntry {
        five_prime_pos: u32,
        idx: usize,
        read_support: usize,
    }

    #[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
    struct SlSupportKey<'a> {
        chrom: &'a str,
        strand: Strand,
        junctions: &'a [u32],
    }

    let mut groups: HashMap<SlSupportKey<'_>, Vec<SlSupportEntry>> = HashMap::new();
    for (i, track_i) in tracks.iter().enumerate() {
        if !is_sl_supported_read(track_i, sw_score) {
            continue;
        }
        let Some(five_prime_pos) = five_prime_position(&track_i.tx) else {
            continue;
        };

        groups
            .entry(SlSupportKey {
                chrom: track_i.tx.chrom.as_str(),
                strand: track_i.tx.strand,
                junctions: junctions_cache[i].as_slice(),
            })
            .or_default()
            .push(SlSupportEntry {
                five_prime_pos,
                idx: i,
                read_support: track_read_support(track_i),
            });
    }

    for entries in groups.values_mut() {
        entries.sort_by_key(|entry| (entry.five_prime_pos, entry.idx));

        let mut prefix_support = Vec::with_capacity(entries.len() + 1);
        prefix_support.push(0usize);
        for entry in entries.iter() {
            let previous = *prefix_support.last().expect("prefix has initial zero");
            prefix_support.push(previous + entry.read_support);
        }

        let mut left = 0usize;
        let mut right = 0usize;
        for idx in 0..entries.len() {
            let center = entries[idx].five_prime_pos;
            while entries[left]
                .five_prime_pos
                .saturating_add(sl_options.five_prime_cluster_offset)
                < center
            {
                left += 1;
            }
            while right < entries.len()
                && entries[right].five_prime_pos
                    <= center.saturating_add(sl_options.five_prime_cluster_offset)
            {
                right += 1;
            }

            support[entries[idx].idx] = prefix_support[right] - prefix_support[left];
        }
    }

    support
}

fn build_three_prime_cluster_support(
    tracks: &[Track],
    junctions_cache: &[Vec<u32>],
    three_prime_options: ThreePrimeMergeOptions,
) -> Vec<usize> {
    let mut support = vec![0usize; tracks.len()];

    #[derive(Clone, Copy, Debug)]
    struct TerminalSupportEntry {
        three_prime_pos: u32,
        idx: usize,
        read_support: usize,
    }

    #[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
    struct TerminalSupportKey<'a> {
        chrom: &'a str,
        strand: Strand,
        junctions: &'a [u32],
    }

    let mut groups: HashMap<TerminalSupportKey<'_>, Vec<TerminalSupportEntry>> = HashMap::new();
    for (idx, track) in tracks.iter().enumerate() {
        if track.is_reference() || junctions_cache[idx].is_empty() {
            continue;
        }
        let Some(three_prime_pos) = three_prime_position(&track.tx) else {
            continue;
        };

        groups
            .entry(TerminalSupportKey {
                chrom: track.tx.chrom.as_str(),
                strand: track.tx.strand,
                junctions: junctions_cache[idx].as_slice(),
            })
            .or_default()
            .push(TerminalSupportEntry {
                three_prime_pos,
                idx,
                read_support: track_read_support(track),
            });
    }

    for entries in groups.values_mut() {
        entries.sort_by_key(|entry| (entry.three_prime_pos, entry.idx));

        let mut prefix_support = Vec::with_capacity(entries.len() + 1);
        prefix_support.push(0usize);
        for entry in entries.iter() {
            let previous = *prefix_support.last().expect("prefix has initial zero");
            prefix_support.push(previous + entry.read_support);
        }

        let mut left = 0usize;
        let mut right = 0usize;
        for idx in 0..entries.len() {
            let center = entries[idx].three_prime_pos;
            while entries[left]
                .three_prime_pos
                .saturating_add(three_prime_options.three_prime_cluster_offset)
                < center
            {
                left += 1;
            }
            while right < entries.len()
                && entries[right].three_prime_pos
                    <= center.saturating_add(three_prime_options.three_prime_cluster_offset)
            {
                right += 1;
            }

            support[entries[idx].idx] = prefix_support[right] - prefix_support[left];
        }
    }

    support
}

fn sl_protected_from_merge(
    short: &Track,
    long: &Track,
    kind: MergeKind,
    sw_score: i64,
    sl_cluster_support: usize,
    sl_options: SlMergeOptions,
) -> bool {
    if !is_sl_supported_read(short, sw_score) {
        return false;
    }
    if sl_cluster_support < sl_options.min_five_prime_cluster_support {
        return false;
    }

    let offset = match kind {
        MergeKind::SameJunction => sl_options.same_junction_five_prime_end_offset,
        MergeKind::FivePrimeTruncation
        | MergeKind::SingleExonContained
        | MergeKind::SingleExonSameFivePrime => sl_options.partial_five_prime_end_offset,
    };

    !five_prime_ends_match(&short.tx, &long.tx, offset)
}

fn three_prime_protected_from_merge(
    short: &Track,
    long: &Track,
    kind: MergeKind,
    three_prime_cluster_support: usize,
    three_prime_options: ThreePrimeMergeOptions,
) -> bool {
    if kind != MergeKind::SameJunction || short.is_reference() {
        return false;
    }
    if three_prime_cluster_support < three_prime_options.min_three_prime_cluster_support {
        return false;
    }

    three_prime_end_delta(&short.tx, &long.tx)
        .is_some_and(|delta| delta > three_prime_options.same_junction_three_prime_end_offset)
}

fn is_single_exon_same_5prime_in(
    single: &Transcript,
    other: &Transcript,
    other_junctions: &[u32],
    sl_options: SlMergeOptions,
) -> bool {
    if !other_junctions.is_empty()
        || !five_prime_ends_match(single, other, sl_options.partial_five_prime_end_offset)
    {
        return false;
    }

    match single.strand {
        Strand::Plus => single.tx_end <= other.tx_end,
        Strand::Minus => single.tx_start >= other.tx_start,
        Strand::Unknown => false,
    }
}

fn single_exon_merge_kind(
    single: &Transcript,
    other: &Transcript,
    other_junctions: &[u32],
    sl_options: SlMergeOptions,
) -> Option<MergeKind> {
    if is_single_exon_in(single, other, other_junctions) {
        return Some(MergeKind::SingleExonContained);
    }
    if is_single_exon_same_5prime_in(single, other, other_junctions, sl_options) {
        return Some(MergeKind::SingleExonSameFivePrime);
    }
    None
}

fn target_is_preferred_container(
    source_idx: usize,
    target_idx: usize,
    kind: MergeKind,
    source_exon_len: u32,
    target_exon_len: u32,
    is_reference: &[bool],
) -> bool {
    if is_reference[source_idx] {
        return false;
    }
    if is_reference[target_idx] {
        return true;
    }

    match kind {
        MergeKind::SameJunction => {
            target_exon_len > source_exon_len
                || (target_exon_len == source_exon_len && target_idx > source_idx)
        }
        MergeKind::FivePrimeTruncation
        | MergeKind::SingleExonContained
        | MergeKind::SingleExonSameFivePrime => true,
    }
}

struct MergeContext<'a> {
    tracks: &'a [Track],
    junctions_cache: &'a [Vec<u32>],
    exon_lens: &'a [u32],
    is_reference: &'a [bool],
    sl_cluster_support: &'a [usize],
    three_prime_cluster_support: &'a [usize],
    sw_score: i64,
    sl_options: SlMergeOptions,
    three_prime_options: ThreePrimeMergeOptions,
    same_junction_offset: u32,
}

impl MergeContext<'_> {
    fn should_merge_into(&self, source_idx: usize, target_idx: usize) -> bool {
        if source_idx == target_idx {
            return false;
        }

        let source_junctions = &self.junctions_cache[source_idx];
        let kind = if source_junctions.is_empty() {
            let target_junctions = &self.junctions_cache[target_idx];
            let Some(kind) = single_exon_merge_kind(
                &self.tracks[source_idx].tx,
                &self.tracks[target_idx].tx,
                target_junctions,
                self.sl_options,
            ) else {
                return false;
            };
            kind
        } else {
            let Some(kind) = junction_merge_kind(
                source_junctions,
                &self.junctions_cache[target_idx],
                self.same_junction_offset,
            ) else {
                return false;
            };
            kind
        };

        target_is_preferred_container(
            source_idx,
            target_idx,
            kind,
            self.exon_lens[source_idx],
            self.exon_lens[target_idx],
            self.is_reference,
        ) && !sl_protected_from_merge(
            &self.tracks[source_idx],
            &self.tracks[target_idx],
            kind,
            self.sw_score,
            self.sl_cluster_support[source_idx],
            self.sl_options,
        ) && !three_prime_protected_from_merge(
            &self.tracks[source_idx],
            &self.tracks[target_idx],
            kind,
            self.three_prime_cluster_support[source_idx],
            self.three_prime_options,
        )
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

fn build_junction_suffix_index<'a>(
    junctions_cache: &'a [Vec<u32>],
    target_eligible: &[bool],
) -> HashMap<&'a [u32], Vec<usize>> {
    let total_suffixes = junctions_cache
        .iter()
        .zip(target_eligible.iter().copied())
        .filter(|(_, eligible)| *eligible)
        .map(|(junctions, _)| junctions.len())
        .sum();
    let mut suffix_index: HashMap<&'a [u32], Vec<usize>> = HashMap::with_capacity(total_suffixes);

    for (idx, junctions) in junctions_cache.iter().enumerate() {
        if !target_eligible[idx] {
            continue;
        }
        for start in 0..junctions.len() {
            suffix_index
                .entry(&junctions[start..])
                .or_default()
                .push(idx);
        }
    }

    suffix_index
}

fn build_junction_length_index(
    junctions_cache: &[Vec<u32>],
    target_eligible: &[bool],
) -> HashMap<usize, Vec<(u32, u32, usize)>> {
    let mut length_index: HashMap<usize, Vec<(u32, u32, usize)>> = HashMap::new();
    for (idx, junctions) in junctions_cache.iter().enumerate() {
        if !target_eligible[idx] || junctions.is_empty() {
            continue;
        }
        let first = junctions[0];
        let last = *junctions.last().expect("non-empty junctions");
        length_index
            .entry(junctions.len())
            .or_default()
            .push((first, last, idx));
    }
    for bucket in length_index.values_mut() {
        bucket.sort_unstable_by_key(|(first, last, idx)| (*first, *last, *idx));
    }
    length_index
}

fn same_length_window_candidates(
    length_index: &HashMap<usize, Vec<(u32, u32, usize)>>,
    junctions: &[u32],
    offset: u32,
) -> Vec<usize> {
    if offset == 0 || junctions.is_empty() {
        return Vec::new();
    }
    let Some(bucket) = length_index.get(&junctions.len()) else {
        return Vec::new();
    };
    let first = junctions[0];
    let last = *junctions.last().expect("non-empty junctions");
    let lo = first.saturating_sub(offset);
    let hi = first.saturating_add(offset);
    let start = bucket.partition_point(|(value, _, _)| *value < lo);
    let mut out = Vec::new();
    for &(value, other_last, idx) in &bucket[start..] {
        if value > hi {
            break;
        }
        if last.abs_diff(other_last) <= offset {
            out.push(idx);
        }
    }
    out
}

struct SingleExonTargetIndex {
    exact_span: HashMap<(u32, u32), Vec<usize>>,
    single_by_start: Vec<(u32, u32, usize)>,
    five_prime: Vec<(u32, usize)>,
    spliced_plus: Vec<(u32, u32, usize)>,
    spliced_minus: Vec<(u32, u32, usize)>,
}

impl SingleExonTargetIndex {
    fn new(tracks: &[Track], junctions_cache: &[Vec<u32>], target_eligible: &[bool]) -> Self {
        let mut exact_span: HashMap<(u32, u32), Vec<usize>> = HashMap::new();
        let mut single_by_start = Vec::new();
        let mut five_prime = Vec::new();
        let mut spliced_plus = Vec::new();
        let mut spliced_minus = Vec::new();

        for (idx, track) in tracks.iter().enumerate() {
            if !target_eligible[idx] {
                continue;
            }
            let start = track.tx.tx_start.get();
            let end = track.tx.tx_end.get();
            if junctions_cache[idx].is_empty() {
                exact_span.entry((start, end)).or_default().push(idx);
                single_by_start.push((start, end, idx));
                if let Some(pos) = five_prime_position(&track.tx) {
                    five_prime.push((pos, idx));
                }
            } else {
                let last_junction = *junctions_cache[idx].last().expect("non-empty");
                match track.tx.strand {
                    Strand::Plus => spliced_plus.push((last_junction, end, idx)),
                    Strand::Minus => spliced_minus.push((start, last_junction, idx)),
                    Strand::Unknown => {}
                }
            }
        }

        single_by_start.sort_unstable_by_key(|(start, end, idx)| (*start, *end, *idx));
        five_prime.sort_unstable_by_key(|(pos, idx)| (*pos, *idx));
        spliced_plus.sort_unstable_by_key(|(last, end, idx)| (*last, *end, *idx));
        spliced_minus.sort_unstable_by_key(|(start, last, idx)| (*start, *last, *idx));

        Self {
            exact_span,
            single_by_start,
            five_prime,
            spliced_plus,
            spliced_minus,
        }
    }

    fn candidates(
        &self,
        source_idx: usize,
        source: &Transcript,
        sl_options: SlMergeOptions,
    ) -> Vec<usize> {
        let start = source.tx_start.get();
        let end = source.tx_end.get();
        let mut seen = HashSet::new();
        let mut out = Vec::new();
        let mut push = |idx: usize| {
            if idx != source_idx && seen.insert(idx) {
                out.push(idx);
            }
        };

        if let Some(ids) = self.exact_span.get(&(start, end)) {
            for &idx in ids {
                push(idx);
            }
        }

        let prefix = self
            .single_by_start
            .partition_point(|(other_start, _, _)| *other_start < start);
        for &(_, other_end, idx) in &self.single_by_start[..prefix] {
            if other_end > end {
                push(idx);
            }
        }

        if let Some(five_prime_pos) = five_prime_position(source) {
            let lo = five_prime_pos.saturating_sub(sl_options.partial_five_prime_end_offset);
            let hi = five_prime_pos.saturating_add(sl_options.partial_five_prime_end_offset);
            let window_start = self.five_prime.partition_point(|(pos, _)| *pos < lo);
            for &(pos, idx) in &self.five_prime[window_start..] {
                if pos > hi {
                    break;
                }
                push(idx);
            }
        }

        match source.strand {
            Strand::Plus => {
                let prefix = self
                    .spliced_plus
                    .partition_point(|(last_junction, _, _)| *last_junction <= start);
                for &(_, tx_end, idx) in &self.spliced_plus[..prefix] {
                    if end <= tx_end {
                        push(idx);
                    }
                }
            }
            Strand::Minus => {
                let prefix = self
                    .spliced_minus
                    .partition_point(|(tx_start, _, _)| *tx_start <= start);
                for &(_, last_junction, idx) in &self.spliced_minus[..prefix] {
                    if end <= last_junction {
                        push(idx);
                    }
                }
            }
            Strand::Unknown => {}
        }

        out
    }
}

fn build_exact_duplicate_representatives<'a>(
    tracks: &'a [Track],
    junctions_cache: &'a [Vec<u32>],
    is_reference: &[bool],
) -> Vec<Option<usize>> {
    #[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
    struct ExactDuplicateKey<'a> {
        chrom: &'a str,
        strand: Strand,
        tx_start: u32,
        tx_end: u32,
        junctions: &'a [u32],
    }

    let mut groups: HashMap<ExactDuplicateKey<'_>, Vec<usize>> = HashMap::new();
    for (idx, track) in tracks.iter().enumerate() {
        if is_reference[idx] {
            continue;
        }

        groups
            .entry(ExactDuplicateKey {
                chrom: track.tx.chrom.as_str(),
                strand: track.tx.strand,
                tx_start: track.tx.tx_start.get(),
                tx_end: track.tx.tx_end.get(),
                junctions: junctions_cache[idx].as_slice(),
            })
            .or_default()
            .push(idx);
    }

    let mut representatives = vec![None; tracks.len()];
    for group in groups.values() {
        if group.len() <= 1 {
            continue;
        }

        // Same-junction equal-length tie breaking keeps the latest read, so use the
        // last exact duplicate as the only non-reference target for this group.
        let representative = *group.last().expect("non-empty group");
        for &idx in group {
            if idx != representative {
                representatives[idx] = Some(representative);
            }
        }
    }

    representatives
}

#[cfg(test)]
fn junction_simple_merge(tracks: &mut [Track], sw_score: i64) -> Vec<usize> {
    junction_simple_merge_with_options(
        tracks,
        sw_score,
        SlMergeOptions::default(),
        ThreePrimeMergeOptions::default(),
        0,
    )
}

fn junction_simple_merge_with_options(
    tracks: &mut [Track],
    sw_score: i64,
    sl_options: SlMergeOptions,
    three_prime_options: ThreePrimeMergeOptions,
    same_junction_offset: u32,
) -> Vec<usize> {
    let junctions_cache: Vec<Vec<u32>> = tracks
        .iter()
        .map(|track| junction_positions(&track.tx))
        .collect();
    let exon_lens: Vec<u32> = tracks.iter().map(|track| exon_len(&track.tx)).collect();
    let is_reference: Vec<bool> = tracks.iter().map(Track::is_reference).collect();
    let exact_duplicate_representative =
        build_exact_duplicate_representatives(tracks, &junctions_cache, &is_reference);
    let target_eligible: Vec<bool> = exact_duplicate_representative
        .iter()
        .map(Option::is_none)
        .collect();
    let sl_cluster_support =
        build_sl_five_prime_cluster_support(tracks, &junctions_cache, sw_score, sl_options);
    let three_prime_cluster_support =
        build_three_prime_cluster_support(tracks, &junctions_cache, three_prime_options);

    let suffix_index = build_junction_suffix_index(&junctions_cache, &target_eligible);
    let length_index = build_junction_length_index(&junctions_cache, &target_eligible);
    let single_exon_index = SingleExonTargetIndex::new(tracks, &junctions_cache, &target_eligible);

    let mut dropped: Vec<bool> = vec![false; tracks.len()];
    for i in 0..tracks.len() {
        if dropped[i] {
            continue;
        }

        let junctions_i = &junctions_cache[i];
        if junctions_i.is_empty() {
            for j in single_exon_index.candidates(i, &tracks[i].tx, sl_options) {
                if dropped[j] && !is_reference[j] {
                    continue;
                }

                let should_merge = {
                    let ctx = MergeContext {
                        tracks,
                        junctions_cache: &junctions_cache,
                        exon_lens: &exon_lens,
                        is_reference: &is_reference,
                        sl_cluster_support: &sl_cluster_support,
                        three_prime_cluster_support: &three_prime_cluster_support,
                        sw_score,
                        sl_options,
                        three_prime_options,
                        same_junction_offset,
                    };
                    ctx.should_merge_into(i, j)
                };

                if should_merge {
                    dropped[i] = true;
                    let (short, long) = get_two_mut(tracks, i, j);
                    long.subreads.extend(short.subreads.iter().cloned());
                }
            }
            continue;
        }

        let exact_candidates = suffix_index.get(junctions_i.as_slice());
        let same_length_candidates =
            same_length_window_candidates(&length_index, junctions_i, same_junction_offset);
        if exact_candidates.is_none() && same_length_candidates.is_empty() {
            continue;
        }

        let mut seen_candidates: HashSet<usize> = HashSet::new();
        for &j in exact_candidates
            .into_iter()
            .flatten()
            .chain(same_length_candidates.iter())
        {
            if !seen_candidates.insert(j) {
                continue;
            }
            if i == j {
                continue;
            }
            if dropped[j] && !is_reference[j] {
                continue;
            }
            debug_assert!(target_eligible[j]);

            let should_merge = {
                let ctx = MergeContext {
                    tracks,
                    junctions_cache: &junctions_cache,
                    exon_lens: &exon_lens,
                    is_reference: &is_reference,
                    sl_cluster_support: &sl_cluster_support,
                    three_prime_cluster_support: &three_prime_cluster_support,
                    sw_score,
                    sl_options,
                    three_prime_options,
                    same_junction_offset,
                };
                ctx.should_merge_into(i, j)
            };

            if should_merge {
                dropped[i] = true;
                let (short, long) = get_two_mut(tracks, i, j);
                long.subreads.extend(short.subreads.iter().cloned());
            }
        }
    }

    let mut keep_vec: Vec<usize> = Vec::with_capacity(tracks.len());
    for (idx, _track) in tracks.iter().enumerate() {
        if !dropped[idx] || is_reference[idx] {
            keep_vec.push(idx);
        }
    }
    keep_vec
}

#[cfg(test)]
fn junction_simple_merge_naive(tracks: &mut [Track], sw_score: i64) -> Vec<usize> {
    junction_simple_merge_naive_with_options(tracks, sw_score, SlMergeOptions::default(), 0)
}

#[cfg(test)]
fn junction_simple_merge_naive_with_options(
    tracks: &mut [Track],
    sw_score: i64,
    sl_options: SlMergeOptions,
    same_junction_offset: u32,
) -> Vec<usize> {
    let junctions_cache: Vec<Vec<u32>> = tracks
        .iter()
        .map(|track| junction_positions(&track.tx))
        .collect();
    let exon_lens: Vec<u32> = tracks.iter().map(|track| exon_len(&track.tx)).collect();
    let is_reference: Vec<bool> = tracks.iter().map(Track::is_reference).collect();
    let sl_cluster_support =
        build_sl_five_prime_cluster_support(tracks, &junctions_cache, sw_score, sl_options);
    let three_prime_options = ThreePrimeMergeOptions::default();
    let three_prime_cluster_support =
        build_three_prime_cluster_support(tracks, &junctions_cache, three_prime_options);

    let mut dropped: Vec<bool> = vec![false; tracks.len()];
    for i in 0..tracks.len() {
        if dropped[i] {
            continue;
        }

        for j in 0..tracks.len() {
            if i == j {
                continue;
            }
            if dropped[j] && !is_reference[j] {
                continue;
            }

            let should_merge = {
                let ctx = MergeContext {
                    tracks,
                    junctions_cache: &junctions_cache,
                    exon_lens: &exon_lens,
                    is_reference: &is_reference,
                    sl_cluster_support: &sl_cluster_support,
                    three_prime_cluster_support: &three_prime_cluster_support,
                    sw_score,
                    sl_options,
                    three_prime_options,
                    same_junction_offset,
                };
                ctx.should_merge_into(i, j)
            };

            if should_merge {
                dropped[i] = true;
                let (short, long) = get_two_mut(tracks, i, j);
                long.subreads.extend(short.subreads.iter().cloned());
            }
        }
    }

    let mut keep_vec: Vec<usize> = Vec::with_capacity(tracks.len());
    for (idx, _track) in tracks.iter().enumerate() {
        if !dropped[idx] || is_reference[idx] {
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
    let mut index_by_source_name_and_structure: HashMap<(TrackSource, String, String), usize> =
        HashMap::new();

    for track in tracks {
        // A read label identifies a molecule, not a unique alignment. Preserve
        // structurally distinct alignments and coalesce only exact structural
        // copies of the same source/name.
        let key = (
            track.source,
            track.tx.name.clone(),
            crate::identity::novel_isoform_id(&track.tx),
        );
        if let Some(&idx) = index_by_source_name_and_structure.get(&key) {
            out[idx].subreads.extend(track.subreads);
        } else {
            index_by_source_name_and_structure.insert(key, out.len());
            out.push(track);
        }
    }

    out
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

fn merge_one_batch(
    mut tracks: Vec<Track>,
    sw_score: i64,
    sl_options: SlMergeOptions,
    three_prime_options: ThreePrimeMergeOptions,
    same_junction_offset: u32,
) -> Vec<Track> {
    let keep_indices = junction_simple_merge_with_options(
        &mut tracks,
        sw_score,
        sl_options,
        three_prime_options,
        same_junction_offset,
    );
    select_tracks_by_keep_indices(tracks, keep_indices)
}

fn merge_read_batches(
    ref_tracks: &[Track],
    read_tracks: Vec<Track>,
    batch_size: usize,
    sw_score: i64,
    sl_options: SlMergeOptions,
    three_prime_options: ThreePrimeMergeOptions,
    same_junction_offset: u32,
) -> (Vec<Track>, bool) {
    let mut anchors: Vec<Track> = ref_tracks.to_vec();
    let mut changed = false;

    for chunk in read_tracks.chunks(batch_size.max(1)) {
        let mut batch: Vec<Track> = Vec::with_capacity(anchors.len() + chunk.len());
        batch.extend(anchors);
        batch.extend(chunk.iter().cloned());

        let before_len = batch.len();
        let merged = merge_one_batch(
            batch,
            sw_score,
            sl_options,
            three_prime_options,
            same_junction_offset,
        );
        if merged.len() < before_len {
            changed = true;
        }

        anchors = merge_tracks_by_name(merged);
    }

    (anchors, changed)
}

fn batch_junction_simple_merge(
    tracks: Vec<Track>,
    sw_score: i64,
    batch_size: usize,
    max_rounds: usize,
    sl_options: SlMergeOptions,
    three_prime_options: ThreePrimeMergeOptions,
    same_junction_offset: u32,
) -> Vec<Track> {
    let batch_size = batch_size.max(1);
    let max_rounds = max_rounds.max(1);

    let mut tracks = tracks;
    let mut rounds = 0usize;
    let mut previous_len = tracks.len();

    // Read batching always keeps references available as potential containers.
    // Use the caller's SW cutoff for every round so negative cutoffs disable only
    // SW/SL 5' protection, not ordinary truncation merging.
    while rounds < max_rounds {
        let (refs, reads) = split_reference_and_read_tracks(tracks);
        if reads.len() <= batch_size {
            let mut combined: Vec<Track> = Vec::with_capacity(refs.len() + reads.len());
            combined.extend(refs);
            combined.extend(reads);
            return merge_one_batch(
                combined,
                sw_score,
                sl_options,
                three_prime_options,
                same_junction_offset,
            );
        }

        let (merged, changed) = merge_read_batches(
            &refs,
            reads,
            batch_size,
            sw_score,
            sl_options,
            three_prime_options,
            same_junction_offset,
        );
        tracks = merged;

        if !changed || tracks.len() >= previous_len {
            break;
        }
        previous_len = tracks.len();
        rounds += 1;
    }

    let (refs, reads) = split_reference_and_read_tracks(tracks);
    let mut combined: Vec<Track> = Vec::with_capacity(refs.len() + reads.len());
    combined.extend(refs);
    combined.extend(reads);
    merge_one_batch(
        combined,
        sw_score,
        sl_options,
        three_prime_options,
        same_junction_offset,
    )
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

fn split_tracks_into_loci(tracks: Vec<Track>) -> Vec<Vec<Track>> {
    if tracks.is_empty() {
        return Vec::new();
    }
    if tracks.len() == 1 {
        return vec![tracks];
    }

    let mut tracks: Vec<Option<Track>> = tracks.into_iter().map(Some).collect();
    let mut order: Vec<usize> = (0..tracks.len()).collect();
    order.sort_by(|left, right| {
        let left_tx = &tracks[*left].as_ref().expect("track present").tx;
        let right_tx = &tracks[*right].as_ref().expect("track present").tx;
        left_tx
            .chrom
            .cmp(&right_tx.chrom)
            .then_with(|| left_tx.tx_start.cmp(&right_tx.tx_start))
            .then_with(|| left_tx.tx_end.cmp(&right_tx.tx_end))
            .then_with(|| left_tx.strand.cmp(&right_tx.strand))
    });

    let sorted_records: Vec<Transcript> = order
        .iter()
        .map(|idx| tracks[*idx].as_ref().expect("track present").tx.clone())
        .collect();
    let loci = cluster_by_span(&sorted_records, StrandMode::Match);

    let mut out: Vec<Vec<Track>> = Vec::with_capacity(loci.len());
    for locus in loci {
        let mut locus_tracks: Vec<Track> = Vec::with_capacity(locus.members.len());
        for member in locus.members {
            let original_idx = order[member];
            locus_tracks.push(
                tracks[original_idx]
                    .take()
                    .expect("track already consumed for locus"),
            );
        }
        out.push(locus_tracks);
    }

    out
}

struct PartitionResult {
    isoforms: Vec<Transcript>,
    pairs: Vec<(String, String)>,
    represented_read_indices: HashSet<usize>,
    rare_read_indices: Vec<usize>,
    unmatched_read_indices: Vec<usize>,
    downsampled_read_indices: Vec<usize>,
}

struct WorkItem {
    index: usize,
    key: PartitionKey,
    ref_indices: Vec<usize>,
    read_indices: Vec<usize>,
}

#[derive(Default)]
struct PartitionWorkerState {
    key: Option<String>,
    started_at: Option<Instant>,
}

fn seed_for_locus(base_seed: u64, tracks: &[Track]) -> u64 {
    let mut hash = crate::rng::fnv1a64(b"clusterj-locus");
    if let Some(first) = tracks.first() {
        crate::rng::update_fnv1a64(&mut hash, first.tx.chrom.as_bytes());
        crate::rng::update_fnv1a64(&mut hash, &[first.tx.strand.as_char() as u8]);
    }
    let start = tracks
        .iter()
        .map(|track| track.tx.tx_start.get())
        .min()
        .unwrap_or(0);
    let end = tracks
        .iter()
        .map(|track| track.tx.tx_end.get())
        .max()
        .unwrap_or(0);
    crate::rng::update_fnv1a64(&mut hash, &start.to_le_bytes());
    crate::rng::update_fnv1a64(&mut hash, &end.to_le_bytes());
    crate::rng::update_fnv1a64(&mut hash, &(tracks.len() as u64).to_le_bytes());
    base_seed ^ hash
}

fn downsample_locus_tracks(
    tracks: Vec<Track>,
    runtime: ClusterjRuntimeOptions,
    downsampled_read_indices: &mut Vec<usize>,
) -> Vec<Track> {
    if runtime.max_reads_per_locus == 0 {
        return tracks;
    }

    let (refs, reads) = split_reference_and_read_tracks(tracks);
    let original_reads = reads.len();
    if original_reads <= runtime.max_reads_per_locus {
        let mut combined = refs;
        combined.extend(reads);
        return combined;
    }

    let seed = seed_for_locus(runtime.downsample_seed, &reads);
    let keep =
        crate::rng::reservoir_sample_indices(original_reads, runtime.max_reads_per_locus, seed);
    let mut keep_mask = vec![false; original_reads];
    for idx in &keep {
        keep_mask[*idx] = true;
    }

    let mut sampled = Vec::with_capacity(keep.len());
    for (idx, track) in reads.into_iter().enumerate() {
        if keep_mask[idx] {
            sampled.push(track);
        } else {
            downsampled_read_indices.extend(track.subreads.iter().map(|subread| subread.index));
        }
    }

    eprintln!(
        "clusterj: subsample locus chrom={} strand={} original_reads={} sampled_reads={} seed={}",
        refs.first()
            .or(sampled.first())
            .map(|track| track.tx.chrom.as_str())
            .unwrap_or("none"),
        refs.first()
            .or(sampled.first())
            .map(|track| track.tx.strand.as_char())
            .unwrap_or('.'),
        original_reads,
        sampled.len(),
        seed
    );

    let mut combined = refs;
    combined.extend(sampled);
    combined
}

#[allow(clippy::too_many_arguments)]
fn process_partition(
    references: &[Transcript],
    reads: &[Transcript],
    ref_indices: &[usize],
    read_indices: &[usize],
    sw_score: i64,
    batch_size: usize,
    batch_rounds: usize,
    name2_mode: Name2Mode,
    sl_options: SlMergeOptions,
    three_prime_options: ThreePrimeMergeOptions,
    junction_correction: JunctionCorrectionOptions,
    runtime: ClusterjRuntimeOptions,
) -> PartitionResult {
    let mut tracks: Vec<Track> = Vec::with_capacity(ref_indices.len() + read_indices.len());
    for &idx in ref_indices {
        tracks.push(Track::reference(references[idx].clone()));
    }
    for &idx in read_indices {
        tracks.push(Track::read(reads[idx].clone(), idx));
    }

    let mut kept: Vec<Track> = Vec::new();
    let mut rare_read_indices: Vec<usize> = Vec::new();
    let mut unmatched_read_indices: Vec<usize> = Vec::new();
    let mut downsampled_read_indices: Vec<usize> = Vec::new();

    // Split before junction correction so a read at a disjoint locus cannot borrow junction
    // support from an unrelated reference on the same chromosome and strand. This matches the
    // overlap-mode contract: loci with no reference anchor are returned as unused.
    for locus_tracks in split_tracks_into_loci(tracks) {
        if !locus_tracks.iter().any(Track::is_reference) {
            unmatched_read_indices.extend(
                locus_tracks
                    .iter()
                    .flat_map(|track| track.subreads.iter().map(|subread| subread.index)),
            );
            continue;
        }

        let locus_tracks =
            downsample_locus_tracks(locus_tracks, runtime, &mut downsampled_read_indices);
        let (corrected, rare) = flow_junction_correct(
            locus_tracks,
            junction_correction.min_support,
            junction_correction.offset,
        );
        rare_read_indices.extend(
            rare.iter()
                .flat_map(|track| track.subreads.iter().map(|subread| subread.index)),
        );

        for mut corrected_locus in split_tracks_into_loci(corrected) {
            let mut locus_kept = if batch_size == 0 {
                let keep_indices = junction_simple_merge_with_options(
                    &mut corrected_locus,
                    sw_score,
                    sl_options,
                    three_prime_options,
                    junction_correction.offset,
                );
                select_tracks_by_keep_indices(corrected_locus, keep_indices)
            } else {
                batch_junction_simple_merge(
                    corrected_locus,
                    sw_score,
                    batch_size,
                    batch_rounds,
                    sl_options,
                    three_prime_options,
                    junction_correction.offset,
                )
            };
            kept.append(&mut locus_kept);
        }
    }

    let represented_read_indices: HashSet<usize> = kept
        .iter()
        .flat_map(|track| track.subreads.iter().map(|subread| subread.index))
        .collect();
    for track in &mut kept {
        if track.is_read() {
            track.tx.name = crate::identity::novel_isoform_id(&track.tx);
        }
    }
    // Two independently retained representatives with the same gene and exact
    // structure are one catalog isoform. Structural IDs make that equality
    // explicit, so coalesce their molecule memberships before serialization.
    let mut kept = merge_tracks_by_name(kept);
    kept.sort_by(|left, right| left.tx.name.cmp(&right.tx.name));
    update_name2(&mut kept, name2_mode);
    let pairs = build_read_to_isoform(&kept);
    let isoforms = kept.into_iter().map(|track| track.tx).collect();

    PartitionResult {
        isoforms,
        pairs,
        represented_read_indices,
        rare_read_indices,
        unmatched_read_indices,
        downsampled_read_indices,
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
    try_clusterj(
        reads,
        references,
        threads,
        sw_score,
        batch_size,
        batch_rounds,
    )
    .unwrap_or_else(|error| panic!("invalid junction-clustering options: {error}"))
}

/// Junction-cluster with default options, returning invalid configuration errors.
pub fn try_clusterj(
    reads: &[Transcript],
    references: Option<&[Transcript]>,
    threads: usize,
    sw_score: i64,
    batch_size: usize,
    batch_rounds: usize,
) -> Result<ClusterResult, crate::config::ParameterError> {
    try_clusterj_with_name2_mode(
        reads,
        references,
        threads,
        sw_score,
        batch_size,
        batch_rounds,
        Name2Mode::Full,
    )
}

pub fn clusterj_with_name2_mode(
    reads: &[Transcript],
    references: Option<&[Transcript]>,
    threads: usize,
    sw_score: i64,
    batch_size: usize,
    batch_rounds: usize,
    name2_mode: Name2Mode,
) -> ClusterResult {
    try_clusterj_with_name2_mode(
        reads,
        references,
        threads,
        sw_score,
        batch_size,
        batch_rounds,
        name2_mode,
    )
    .unwrap_or_else(|error| panic!("invalid junction-clustering options: {error}"))
}

/// Junction-cluster with an explicit name2 mode, returning invalid configuration errors.
#[allow(clippy::too_many_arguments)]
pub fn try_clusterj_with_name2_mode(
    reads: &[Transcript],
    references: Option<&[Transcript]>,
    threads: usize,
    sw_score: i64,
    batch_size: usize,
    batch_rounds: usize,
    name2_mode: Name2Mode,
) -> Result<ClusterResult, crate::config::ParameterError> {
    try_clusterj_with_options(
        reads,
        references,
        threads,
        sw_score,
        batch_size,
        batch_rounds,
        name2_mode,
        SlMergeOptions::default(),
        ThreePrimeMergeOptions::default(),
        JunctionCorrectionOptions::default(),
    )
}

#[allow(clippy::too_many_arguments)]
pub fn clusterj_with_options(
    reads: &[Transcript],
    references: Option<&[Transcript]>,
    threads: usize,
    sw_score: i64,
    batch_size: usize,
    batch_rounds: usize,
    name2_mode: Name2Mode,
    sl_options: SlMergeOptions,
    three_prime_options: ThreePrimeMergeOptions,
    junction_correction: JunctionCorrectionOptions,
) -> ClusterResult {
    try_clusterj_with_options(
        reads,
        references,
        threads,
        sw_score,
        batch_size,
        batch_rounds,
        name2_mode,
        sl_options,
        three_prime_options,
        junction_correction,
    )
    .unwrap_or_else(|error| panic!("invalid junction-clustering options: {error}"))
}

/// Junction-cluster with explicit options, returning invalid configuration errors.
#[allow(clippy::too_many_arguments)]
pub fn try_clusterj_with_options(
    reads: &[Transcript],
    references: Option<&[Transcript]>,
    threads: usize,
    sw_score: i64,
    batch_size: usize,
    batch_rounds: usize,
    name2_mode: Name2Mode,
    sl_options: SlMergeOptions,
    three_prime_options: ThreePrimeMergeOptions,
    junction_correction: JunctionCorrectionOptions,
) -> Result<ClusterResult, crate::config::ParameterError> {
    let (result, summary) = try_clusterj_with_options_and_summary(
        reads,
        references,
        threads,
        sw_score,
        batch_size,
        batch_rounds,
        name2_mode,
        sl_options,
        three_prime_options,
        junction_correction,
    )?;
    summary.emit();
    Ok(result)
}

#[allow(clippy::too_many_arguments)]
pub fn clusterj_with_options_and_summary(
    reads: &[Transcript],
    references: Option<&[Transcript]>,
    threads: usize,
    sw_score: i64,
    batch_size: usize,
    batch_rounds: usize,
    name2_mode: Name2Mode,
    sl_options: SlMergeOptions,
    three_prime_options: ThreePrimeMergeOptions,
    junction_correction: JunctionCorrectionOptions,
) -> (ClusterResult, JunctionClusterSummary) {
    try_clusterj_with_options_and_summary(
        reads,
        references,
        threads,
        sw_score,
        batch_size,
        batch_rounds,
        name2_mode,
        sl_options,
        three_prime_options,
        junction_correction,
    )
    .unwrap_or_else(|error| panic!("invalid junction-clustering options: {error}"))
}

/// Junction-cluster with summary output, returning invalid configuration errors.
#[allow(clippy::too_many_arguments)]
pub fn try_clusterj_with_options_and_summary(
    reads: &[Transcript],
    references: Option<&[Transcript]>,
    threads: usize,
    sw_score: i64,
    batch_size: usize,
    batch_rounds: usize,
    name2_mode: Name2Mode,
    sl_options: SlMergeOptions,
    three_prime_options: ThreePrimeMergeOptions,
    junction_correction: JunctionCorrectionOptions,
) -> Result<(ClusterResult, JunctionClusterSummary), crate::config::ParameterError> {
    try_clusterj_with_runtime_options_and_summary(
        reads,
        references,
        threads,
        sw_score,
        batch_size,
        batch_rounds,
        name2_mode,
        sl_options,
        three_prime_options,
        junction_correction,
        ClusterjRuntimeOptions::default(),
    )
}

/// Junction-cluster with summary output and explicit runtime bounds.
#[allow(clippy::too_many_arguments)]
pub fn try_clusterj_with_runtime_options_and_summary(
    reads: &[Transcript],
    references: Option<&[Transcript]>,
    threads: usize,
    sw_score: i64,
    batch_size: usize,
    batch_rounds: usize,
    name2_mode: Name2Mode,
    sl_options: SlMergeOptions,
    three_prime_options: ThreePrimeMergeOptions,
    junction_correction: JunctionCorrectionOptions,
    runtime: ClusterjRuntimeOptions,
) -> Result<(ClusterResult, JunctionClusterSummary), crate::config::ParameterError> {
    let threads = crate::config::WorkerThreads::new(threads)?.get();
    crate::config::BatchRounds::new(batch_rounds)?;
    sl_options.validate()?;
    three_prime_options.validate()?;
    junction_correction.validate()?;
    crate::identity::validate_read_ids(reads)
        .map_err(crate::config::ParameterError::invalid_identity)?;
    let references = match references {
        Some(references) => references,
        None => {
            let summary = JunctionClusterSummary {
                input_reads: reads.len(),
                represented_reads: 0,
                mapping_rows: 0,
                rare_reads: 0,
                unmatched_reads: reads.len(),
                unused_reads: reads.len(),
            };
            return Ok((
                ClusterResult {
                    isoforms: Vec::new(),
                    read_to_isoform: Vec::new(),
                    unused: reads.to_vec(),
                },
                summary,
            ));
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
    let mut unmatched_read_indices: Vec<usize> = Vec::new();
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
    let mut represented_read_indices: HashSet<usize> = HashSet::new();
    let mut rare_read_indices: Vec<usize> = Vec::new();

    let mut keys: Vec<PartitionKey> = refs_by_key.keys().cloned().collect();
    keys.sort_by(|a, b| a.chrom.cmp(&b.chrom).then_with(|| a.strand.cmp(&b.strand)));

    let mut work: Vec<WorkItem> = Vec::with_capacity(keys.len());
    for (index, key) in keys.iter().enumerate() {
        work.push(WorkItem {
            index,
            key: key.clone(),
            ref_indices: refs_by_key.remove(key).unwrap_or_default(),
            read_indices: reads_by_key.remove(key).unwrap_or_default(),
        });
    }

    let started = Instant::now();
    let total = keys.len();
    let done = Arc::new(AtomicUsize::new(0));
    let worker_count = if threads == 1 || work.len() <= 1 {
        1
    } else {
        threads.min(keys.len().max(1))
    };
    let worker_states: Arc<Vec<Mutex<PartitionWorkerState>>> = Arc::new(
        (0..worker_count)
            .map(|_| Mutex::new(PartitionWorkerState::default()))
            .collect(),
    );
    let (heartbeat_stop_tx, heartbeat_handle) = if runtime.heartbeat_seconds > 0 && total > 0 {
        use std::sync::mpsc::RecvTimeoutError;

        let (stop_tx, stop_rx) = std::sync::mpsc::channel::<()>();
        let done = Arc::clone(&done);
        let worker_states = Arc::clone(&worker_states);
        let heartbeat_seconds = runtime.heartbeat_seconds;
        let heartbeat_top = runtime.heartbeat_top.max(1);
        let handle = std::thread::spawn(move || {
            let mut last_done = done.load(Ordering::Relaxed);
            loop {
                match stop_rx.recv_timeout(Duration::from_secs(heartbeat_seconds)) {
                    Ok(()) | Err(RecvTimeoutError::Disconnected) => break,
                    Err(RecvTimeoutError::Timeout) => {}
                }

                let done_now = done.load(Ordering::Relaxed);
                eprintln!(
                    "heartbeat: {done_now}/{total} elapsed={:?}",
                    started.elapsed()
                );
                if done_now == last_done && done_now < total {
                    let mut inflight = Vec::new();
                    for state_lock in worker_states.iter() {
                        let Ok(state) = state_lock.lock() else {
                            continue;
                        };
                        let (Some(key), Some(started_at)) = (state.key.as_ref(), state.started_at)
                        else {
                            continue;
                        };
                        inflight.push((key.clone(), started_at.elapsed()));
                    }
                    if inflight.is_empty() {
                        eprintln!("heartbeat: no in-flight partitions (all workers idle?)");
                    } else {
                        inflight.sort_by(|left, right| right.1.cmp(&left.1));
                        let mut line = String::from("in_flight(top):");
                        for (key, duration) in inflight.into_iter().take(heartbeat_top) {
                            line.push(' ');
                            line.push_str(&format!("{key}={:.1}s", duration.as_secs_f64()));
                        }
                        eprintln!("{line}");
                    }
                }
                last_done = done_now;
                if done_now >= total {
                    break;
                }
            }
        });
        (Some(stop_tx), Some(handle))
    } else {
        (None, None)
    };

    let mut parts: Vec<Option<PartitionResult>> = (0..keys.len()).map(|_| None).collect();
    if threads == 1 || work.len() <= 1 {
        for item in work {
            let label = format!("{}:{}", item.key.chrom, item.key.strand.as_char());
            if let Ok(mut state) = worker_states[0].lock() {
                *state = PartitionWorkerState {
                    key: Some(label),
                    started_at: Some(Instant::now()),
                };
            }
            parts[item.index] = Some(process_partition(
                references,
                reads,
                &item.ref_indices,
                &item.read_indices,
                sw_score,
                batch_size,
                batch_rounds,
                name2_mode,
                sl_options,
                three_prime_options,
                junction_correction,
                runtime,
            ));
            if let Ok(mut state) = worker_states[0].lock() {
                *state = PartitionWorkerState::default();
            }
            done.fetch_add(1, Ordering::Relaxed);
        }
    } else {
        let queue = Arc::new(Mutex::new(work));
        let (tx, rx) = std::sync::mpsc::channel::<(usize, PartitionResult)>();

        std::thread::scope(|scope| {
            for worker_idx in 0..worker_count {
                let queue = Arc::clone(&queue);
                let tx = tx.clone();
                let done = Arc::clone(&done);
                let worker_states = Arc::clone(&worker_states);

                scope.spawn(move || loop {
                    let item = {
                        let mut guard = queue.lock().expect("work queue poisoned");
                        guard.pop()
                    };
                    let Some(item) = item else {
                        break;
                    };

                    let label = format!("{}:{}", item.key.chrom, item.key.strand.as_char());
                    if let Ok(mut state) = worker_states[worker_idx].lock() {
                        *state = PartitionWorkerState {
                            key: Some(label),
                            started_at: Some(Instant::now()),
                        };
                    }
                    let result = process_partition(
                        references,
                        reads,
                        &item.ref_indices,
                        &item.read_indices,
                        sw_score,
                        batch_size,
                        batch_rounds,
                        name2_mode,
                        sl_options,
                        three_prime_options,
                        junction_correction,
                        runtime,
                    );
                    if let Ok(mut state) = worker_states[worker_idx].lock() {
                        *state = PartitionWorkerState::default();
                    }
                    done.fetch_add(1, Ordering::Relaxed);
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

    if let Some(tx) = heartbeat_stop_tx {
        let _ = tx.send(());
    }
    if let Some(handle) = heartbeat_handle {
        let _ = handle.join();
    }

    let mut downsampled_read_indices: Vec<usize> = Vec::new();
    for part in parts.into_iter().flatten() {
        all_isoforms.extend(part.isoforms);
        all_pairs.extend(part.pairs);
        represented_read_indices.extend(part.represented_read_indices);
        rare_read_indices.extend(part.rare_read_indices);
        unmatched_read_indices.extend(part.unmatched_read_indices);
        downsampled_read_indices.extend(part.downsampled_read_indices);
    }

    rare_read_indices.sort_unstable();
    unmatched_read_indices.sort_unstable();
    downsampled_read_indices.sort_unstable();
    let mut unused_read_indices = rare_read_indices.clone();
    unused_read_indices.extend(unmatched_read_indices.iter().copied());
    unused_read_indices.extend(downsampled_read_indices.iter().copied());
    unused_read_indices.sort_unstable();

    let unused_instance_set: HashSet<usize> = unused_read_indices.iter().copied().collect();
    assert_eq!(
        unused_instance_set.len(),
        unused_read_indices.len(),
        "junction clustering classified a read instance as unused more than once"
    );
    assert!(
        represented_read_indices.is_disjoint(&unused_instance_set),
        "junction clustering classified a read instance as both represented and unused"
    );
    assert!(
        represented_read_indices
            .iter()
            .chain(unused_instance_set.iter())
            .all(|idx| *idx < reads.len()),
        "junction clustering produced an out-of-range read instance"
    );
    assert_eq!(
        represented_read_indices.len() + unused_instance_set.len(),
        reads.len(),
        "junction clustering violated read-instance conservation"
    );

    all_pairs.sort_by(|left, right| left.0.cmp(&right.0).then_with(|| left.1.cmp(&right.1)));
    all_isoforms.sort_by(|left, right| left.name.cmp(&right.name));
    crate::identity::validate_isoform_ids(&all_isoforms)
        .map_err(crate::config::ParameterError::invalid_identity)?;
    let mut all_unused: Vec<Transcript> = unused_read_indices
        .iter()
        .map(|&idx| reads[idx].clone())
        .collect();
    all_unused.sort_by(crate::identity::transcript_order);
    let summary = JunctionClusterSummary {
        input_reads: reads.len(),
        represented_reads: represented_read_indices.len(),
        mapping_rows: all_pairs.len(),
        rare_reads: rare_read_indices.len(),
        unmatched_reads: unmatched_read_indices.len(),
        unused_reads: unused_read_indices.len(),
    };

    Ok((
        ClusterResult {
            isoforms: all_isoforms,
            read_to_isoform: all_pairs,
            unused: all_unused,
        },
        summary,
    ))
}

#[cfg(test)]
mod tests {
    use std::sync::atomic::{AtomicUsize, Ordering};

    use proptest::prelude::*;

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
                let instance = ReadInstance {
                    index: read_index,
                    name: format!("read_{read_index:02}"),
                };
                for track in tracks.iter_mut().take(read_index + 1) {
                    track.subreads.insert(instance.clone());
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

    fn make_plain_tx_on(
        chrom: &str,
        name: &str,
        strand: Strand,
        exons: &[(u32, u32)],
        score: u32,
    ) -> Transcript {
        let tx_start = exons.iter().map(|(start, _)| *start).min().unwrap_or(0);
        let tx_end = exons.iter().map(|(_, end)| *end).max().unwrap_or(0);
        let exons = exons
            .iter()
            .map(|(start, end)| Interval::new(Coord::new(*start), Coord::new(*end)).unwrap())
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
                thick_start: Coord::new(tx_start),
                thick_end: Coord::new(tx_end),
                item_rgb: "0".to_owned(),
                extra_fields: Vec::new(),
            },
        )
        .unwrap()
    }

    fn make_track(
        name: &str,
        strand: Strand,
        exons: &[(u32, u32)],
        ttype: &str,
        score: u32,
    ) -> Track {
        static NEXT_READ_INSTANCE: AtomicUsize = AtomicUsize::new(0);

        let tx = make_tx(name, strand, exons, ttype, score);
        if ttype == "isoform_anno" {
            Track::reference(tx)
        } else {
            Track::read(tx, NEXT_READ_INSTANCE.fetch_add(1, Ordering::Relaxed))
        }
    }

    fn has_subread(track: &Track, name: &str) -> bool {
        track.subreads.iter().any(|subread| subread.name == name)
    }

    fn cluster_with_summary(
        reads: &[Transcript],
        references: Option<&[Transcript]>,
    ) -> (ClusterResult, JunctionClusterSummary) {
        clusterj_with_options_and_summary(
            reads,
            references,
            1,
            DEFAULT_SW_SCORE,
            0,
            1,
            Name2Mode::Full,
            SlMergeOptions::default(),
            ThreePrimeMergeOptions::default(),
            JunctionCorrectionOptions::default(),
        )
    }

    fn mapped_isoform_id<'a>(result: &'a ClusterResult, read_id: &str) -> &'a str {
        result
            .read_to_isoform
            .iter()
            .find_map(|(mapped_read, isoform_id)| {
                (mapped_read == read_id).then_some(isoform_id.as_str())
            })
            .unwrap_or_else(|| panic!("missing mapping for read {read_id:?}"))
    }

    fn decoded_subreads(tx: &Transcript) -> HashSet<String> {
        crate::identity::decode_name2(
            tx.extra_fields
                .first()
                .map(String::as_str)
                .unwrap_or("none"),
        )
        .unwrap()
        .into_iter()
        .collect()
    }

    #[test]
    fn junction_simple_merge_matches_naive_on_suffix_and_equal_cases() {
        let mut tracks = vec![
            make_track(
                "short_equal",
                Strand::Plus,
                &[(120, 150), (200, 240)],
                "nanopore_read",
                0,
            ),
            make_track(
                "long_equal_1",
                Strand::Plus,
                &[(100, 150), (200, 260)],
                "nanopore_read",
                0,
            ),
            make_track(
                "long_equal_2",
                Strand::Plus,
                &[(90, 150), (200, 270)],
                "nanopore_read",
                0,
            ),
            make_track(
                "short_suffix",
                Strand::Plus,
                &[(200, 250), (300, 350)],
                "nanopore_read",
                0,
            ),
            make_track(
                "long_suffix",
                Strand::Plus,
                &[(100, 150), (200, 250), (300, 350)],
                "nanopore_read",
                0,
            ),
        ];

        let mut naive_tracks = tracks.clone();
        let keep_naive = junction_simple_merge_naive(&mut naive_tracks, 11);

        let keep_indexed = junction_simple_merge(&mut tracks, 11);

        assert_eq!(keep_indexed, keep_naive);
        assert_eq!(tracks.len(), naive_tracks.len());
        for (indexed, naive) in tracks.iter().zip(naive_tracks.iter()) {
            assert_eq!(indexed.subreads, naive.subreads);
        }
    }

    #[test]
    fn exact_duplicate_target_pruning_matches_naive_kept_outputs() {
        let mut tracks = vec![
            make_track(
                "dup_a",
                Strand::Plus,
                &[(100, 110), (120, 130), (140, 150)],
                "nanopore_read",
                0,
            ),
            make_track(
                "longer_same_junction",
                Strand::Plus,
                &[(90, 110), (120, 130), (140, 160)],
                "nanopore_read",
                0,
            ),
            make_track(
                "dup_b",
                Strand::Plus,
                &[(100, 110), (120, 130), (140, 150)],
                "nanopore_read",
                0,
            ),
            make_track(
                "ref_same_junction",
                Strand::Plus,
                &[(95, 110), (120, 130), (140, 155)],
                "isoform_anno",
                100,
            ),
        ];

        let mut naive_tracks = tracks.clone();
        let keep_naive = junction_simple_merge_naive(&mut naive_tracks, 11);
        let keep_indexed = junction_simple_merge(&mut tracks, 11);

        assert_eq!(keep_indexed, keep_naive);
        for idx in keep_indexed {
            assert_eq!(tracks[idx].subreads, naive_tracks[idx].subreads);
        }
    }

    #[test]
    fn fuzzy_same_junction_merge_uses_offset() {
        let tracks = vec![
            make_track(
                "long_chain",
                Strand::Plus,
                &[(100, 150), (200, 250), (300, 350)],
                "nanopore_read",
                0,
            ),
            make_track(
                "near_chain",
                Strand::Plus,
                &[(102, 154), (204, 249), (304, 348)],
                "nanopore_read",
                0,
            ),
        ];

        let mut exact_tracks = tracks.clone();
        let keep_exact = junction_simple_merge_with_options(
            &mut exact_tracks,
            11,
            SlMergeOptions::default(),
            ThreePrimeMergeOptions::default(),
            0,
        );
        assert_eq!(keep_exact, vec![0, 1]);

        let mut fuzzy_tracks = tracks;
        let keep_fuzzy = junction_simple_merge_with_options(
            &mut fuzzy_tracks,
            11,
            SlMergeOptions::default(),
            ThreePrimeMergeOptions::default(),
            5,
        );
        assert_eq!(keep_fuzzy, vec![0]);
        assert!(has_subread(&fuzzy_tracks[0], "long_chain"));
        assert!(has_subread(&fuzzy_tracks[0], "near_chain"));
    }

    #[test]
    fn fuzzy_boundary_comparison_is_one_to_one() {
        let query = [99, 101, 111];
        let reference = [100, 110, 112];

        assert!(!junctions_equal(&query, &reference, 2));
        let (missed, extra) = compare_ei_by_boundary(&query, &reference, 2);
        assert_eq!(missed.len(), 1);
        assert_eq!(extra.len(), 1);
    }

    #[test]
    fn clusterj_fuzzy_same_junction_uses_junction_correction_offset() {
        let refs = vec![make_tx(
            "ref_anchor",
            Strand::Plus,
            &[(50, 500)],
            "isoform_anno",
            100,
        )];
        let reads = vec![
            make_tx(
                "long_chain",
                Strand::Plus,
                &[(100, 150), (200, 250), (300, 350)],
                "nanopore_read",
                1,
            ),
            make_tx(
                "near_chain",
                Strand::Plus,
                &[(102, 154), (204, 249), (304, 348)],
                "nanopore_read",
                1,
            ),
        ];

        let exact = clusterj_with_options(
            &reads,
            Some(&refs),
            1,
            11,
            0,
            1,
            Name2Mode::Full,
            SlMergeOptions::default(),
            ThreePrimeMergeOptions::default(),
            JunctionCorrectionOptions {
                min_support: 1,
                offset: 0,
            },
        );
        let exact_targets: HashSet<_> = exact
            .read_to_isoform
            .iter()
            .filter(|(read, _)| read == "long_chain" || read == "near_chain")
            .map(|(_, isoform)| isoform.as_str())
            .collect();
        assert_eq!(exact_targets.len(), 2);
        assert!(exact_targets
            .iter()
            .all(|id| id.starts_with(crate::identity::NOVEL_ISOFORM_PREFIX)));

        let fuzzy = clusterj_with_options(
            &reads,
            Some(&refs),
            1,
            11,
            0,
            1,
            Name2Mode::Full,
            SlMergeOptions::default(),
            ThreePrimeMergeOptions::default(),
            JunctionCorrectionOptions {
                min_support: 1,
                offset: 5,
            },
        );
        let fuzzy_targets: HashSet<_> = fuzzy
            .read_to_isoform
            .iter()
            .filter(|(read, _)| read == "long_chain" || read == "near_chain")
            .map(|(_, isoform)| isoform.as_str())
            .collect();
        assert_eq!(fuzzy_targets.len(), 1);
        assert_eq!(
            mapped_isoform_id(&fuzzy, "near_chain"),
            mapped_isoform_id(&fuzzy, "long_chain")
        );
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(64))]

        #[test]
        fn junction_simple_merge_matches_naive_on_random_inputs(
            strand in prop_oneof![Just(Strand::Plus), Just(Strand::Minus)],
            params in prop::collection::vec(
                (
                    1usize..=4,
                    0u32..=500,
                    5u32..=50,
                    0u32..=50,
                    any::<bool>(),
                    0u32..=150,
                ),
                1..=20,
            ),
        ) {
            let mut tracks = Vec::with_capacity(params.len());
            for (idx, (exon_count, tx_start, exon_len, gap_len, is_ref, score)) in params.into_iter().enumerate() {
                let mut exons = Vec::with_capacity(exon_count);
                let mut cursor = tx_start;
                for _ in 0..exon_count {
                    let start = cursor;
                    let end = cursor + exon_len;
                    exons.push((start, end));
                    cursor = end + gap_len;
                }

                let name = format!("t{idx}");
                let ttype = if is_ref { "isoform_anno" } else { "nanopore_read" };
                tracks.push(make_track(&name, strand, &exons, ttype, score));
            }

            let mut naive_tracks = tracks.clone();
            let keep_naive = junction_simple_merge_naive(&mut naive_tracks, 11);

            let keep_indexed = junction_simple_merge(&mut tracks, 11);

            prop_assert_eq!(keep_indexed, keep_naive);
            prop_assert_eq!(tracks.len(), naive_tracks.len());
            for (indexed, naive) in tracks.iter().zip(naive_tracks.iter()) {
                prop_assert_eq!(&indexed.subreads, &naive.subreads);
            }
        }
    }

    #[test]
    fn junction_correction_snaps_low_coverage_sites_and_rebuilds_exons() {
        let reference = make_tx(
            "ref",
            Strand::Plus,
            &[(100, 110), (200, 210)],
            "isoform_anno",
            100,
        );
        let read = make_tx(
            "read",
            Strand::Plus,
            &[(100, 111), (201, 210)],
            "nanopore_read",
            1,
        );

        let tracks = vec![Track::reference(reference), Track::read(read, 0)];
        let (corrected, rare) = flow_junction_correct(tracks, 2, 10);
        assert!(rare.is_empty());

        let corrected_read = corrected
            .iter()
            .find(|track| track.tx.name == "read")
            .expect("corrected read missing");

        let expected = vec![
            Interval::new(Coord::new(100), Coord::new(110)).unwrap(),
            Interval::new(Coord::new(200), Coord::new(210)).unwrap(),
        ];
        assert_eq!(corrected_read.tx.exons, expected);
        assert_eq!(
            corrected_read.tx.metadata().transcript_type(),
            Some("nanopore_read_corrected")
        );
    }

    #[test]
    fn junction_correction_sends_uncorrectable_reads_to_unused() {
        let reference = make_tx(
            "ref",
            Strand::Plus,
            &[(100, 110), (200, 210)],
            "isoform_anno",
            100,
        );
        let read = make_tx(
            "read_uncorrectable",
            Strand::Plus,
            &[(100, 150), (201, 210)],
            "nanopore_read",
            1,
        );

        let tracks = vec![Track::reference(reference), Track::read(read, 0)];
        let (corrected, rare) = flow_junction_correct(tracks, 2, 10);

        assert!(rare
            .iter()
            .any(|track| track.tx.name == "read_uncorrectable"));
        assert!(!corrected
            .iter()
            .any(|track| track.tx.name == "read_uncorrectable"));
    }

    #[test]
    fn junction_correction_rejects_a_snap_outside_the_read_span() {
        let reference = make_tx(
            "ref",
            Strand::Plus,
            &[(80, 90), (200, 220)],
            "isoform_anno",
            100,
        );
        let read = make_tx(
            "read",
            Strand::Plus,
            &[(100, 101), (201, 210)],
            "nanopore_read",
            100,
        );
        let original_exons = read.exons.clone();

        let tracks = vec![Track::reference(reference), Track::read(read, 0)];
        let (corrected, rare) = flow_junction_correct(tracks, 5, 15);

        assert_eq!(
            corrected
                .iter()
                .map(|track| track.tx.name.as_str())
                .collect::<Vec<_>>(),
            ["ref"]
        );
        assert_eq!(rare.len(), 1);
        assert_eq!(rare[0].tx.name, "read");
        assert_eq!(rare[0].tx.exons, original_exons);
    }

    #[test]
    fn corrected_junction_rebuild_rejects_duplicate_boundaries() {
        assert!(rebuild_exons_from_junctions(
            Coord::new(100),
            Coord::new(210),
            &[110, 110, 200, 200],
        )
        .is_none());
    }

    #[test]
    fn clusterj_junction_correction_offset_controls_nearby_snap() {
        let refs = vec![make_tx(
            "ref",
            Strand::Plus,
            &[(100, 110), (200, 250)],
            "isoform_anno",
            100,
        )];
        let reads = vec![make_tx(
            "read_near",
            Strand::Plus,
            &[(100, 123), (213, 250)],
            "nanopore_read",
            1,
        )];

        let default_offset = clusterj_with_options(
            &reads,
            Some(&refs),
            1,
            11,
            0,
            1,
            Name2Mode::Full,
            SlMergeOptions::default(),
            ThreePrimeMergeOptions::default(),
            JunctionCorrectionOptions {
                min_support: 2,
                offset: 10,
            },
        );
        assert!(default_offset
            .unused
            .iter()
            .any(|tx| tx.name == "read_near"));

        let widened_offset = clusterj_with_options(
            &reads,
            Some(&refs),
            1,
            11,
            0,
            1,
            Name2Mode::Full,
            SlMergeOptions::default(),
            ThreePrimeMergeOptions::default(),
            JunctionCorrectionOptions {
                min_support: 2,
                offset: 15,
            },
        );
        assert!(widened_offset.unused.is_empty());
        assert!(widened_offset
            .read_to_isoform
            .iter()
            .any(|(read, _)| read == "read_near"));
    }

    #[test]
    fn clusterj_junction_correction_min_support_controls_novel_junction_retention() {
        let refs = vec![make_tx(
            "ref",
            Strand::Plus,
            &[(100, 110), (200, 250)],
            "isoform_anno",
            100,
        )];
        let reads = vec![
            make_tx(
                "read_novel_a",
                Strand::Plus,
                &[(100, 150), (250, 300)],
                "nanopore_read",
                1,
            ),
            make_tx(
                "read_novel_b",
                Strand::Plus,
                &[(100, 150), (250, 300)],
                "nanopore_read",
                1,
            ),
        ];

        let supported = clusterj_with_options(
            &reads,
            Some(&refs),
            1,
            11,
            0,
            1,
            Name2Mode::Full,
            SlMergeOptions::default(),
            ThreePrimeMergeOptions::default(),
            JunctionCorrectionOptions {
                min_support: 2,
                offset: 10,
            },
        );
        assert!(supported.unused.is_empty());

        let unsupported = clusterj_with_options(
            &reads,
            Some(&refs),
            1,
            11,
            0,
            1,
            Name2Mode::Full,
            SlMergeOptions::default(),
            ThreePrimeMergeOptions::default(),
            JunctionCorrectionOptions {
                min_support: 3,
                offset: 10,
            },
        );
        let unused_names: HashSet<&str> = unsupported
            .unused
            .iter()
            .map(|tx| tx.name.as_str())
            .collect();
        assert!(unused_names.contains("read_novel_a"));
        assert!(unused_names.contains("read_novel_b"));
    }

    #[test]
    fn retained_read_isoform_counts_itself() {
        let refs = vec![make_tx(
            "ref1",
            Strand::Plus,
            &[(100, 200)],
            "isoform_anno",
            100,
        )];
        let reads = vec![make_tx(
            "read1",
            Strand::Plus,
            &[(50, 150)],
            "nanopore_read",
            1,
        )];

        let result = clusterj_with_name2_mode(&reads, Some(&refs), 1, 11, 0, 1, Name2Mode::Full);
        assert!(result.unused.is_empty());

        let stable_id = mapped_isoform_id(&result, "read1").to_owned();
        assert!(stable_id.starts_with(crate::identity::NOVEL_ISOFORM_PREFIX));

        let read_isoform = result
            .isoforms
            .iter()
            .find(|tx| tx.name == stable_id)
            .expect("read isoform missing");
        let name2 = read_isoform
            .extra_fields
            .first()
            .map(String::as_str)
            .unwrap_or("");
        assert!(name2.starts_with("tc_name2_v1:"));
        assert!(decoded_subreads(read_isoform).contains("read1"));

        let counts = crate::count::count_by_subreads(&result.isoforms, &refs).unwrap();
        let read_count = counts
            .iter()
            .find(|record| record.isoform_id == stable_id)
            .expect("missing count record for read isoform");
        assert!((read_count.count - 1.0).abs() < 1e-9);
    }

    #[test]
    fn name2_mode_coverage_writes_only_coverage_but_keeps_mapping() {
        let refs = vec![make_tx(
            "ref1",
            Strand::Plus,
            &[(100, 200)],
            "isoform_anno",
            100,
        )];
        let reads = vec![make_tx(
            "read1",
            Strand::Plus,
            &[(50, 150)],
            "nanopore_read",
            1,
        )];

        let result =
            clusterj_with_name2_mode(&reads, Some(&refs), 1, 11, 0, 1, Name2Mode::Coverage);
        let stable_id = mapped_isoform_id(&result, "read1").to_owned();
        let read_isoform = result
            .isoforms
            .iter()
            .find(|tx| tx.name == stable_id)
            .expect("read isoform missing");
        let name2 = read_isoform
            .extra_fields
            .first()
            .map(String::as_str)
            .unwrap_or("");
        assert!(name2.starts_with('|'));
        assert!(!name2.contains("read1"));
        assert!(stable_id.starts_with(crate::identity::NOVEL_ISOFORM_PREFIX));
    }

    #[test]
    fn name2_mode_none_disables_subread_payload() {
        let refs = vec![make_tx(
            "ref1",
            Strand::Plus,
            &[(100, 200)],
            "isoform_anno",
            100,
        )];
        let reads = vec![make_tx(
            "read1",
            Strand::Plus,
            &[(50, 150)],
            "nanopore_read",
            1,
        )];

        let result = clusterj_with_name2_mode(&reads, Some(&refs), 1, 11, 0, 1, Name2Mode::None);
        let stable_id = mapped_isoform_id(&result, "read1").to_owned();
        let read_isoform = result
            .isoforms
            .iter()
            .find(|tx| tx.name == stable_id)
            .expect("read isoform missing");
        let name2 = read_isoform
            .extra_fields
            .first()
            .map(String::as_str)
            .unwrap_or("");
        assert_eq!(name2, "none");
        assert!(stable_id.starts_with(crate::identity::NOVEL_ISOFORM_PREFIX));
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

        let result = clusterj_with_name2_mode(&reads, Some(&refs), 1, 11, 0, 1, Name2Mode::Full);
        assert_eq!(result.unused.len(), 0);

        let iso_names: HashSet<_> = result.isoforms.iter().map(|tx| tx.name.as_str()).collect();
        assert!(iso_names.contains("ref_a"));
        assert!(iso_names.contains("ref_b"));
        assert!(!iso_names.contains("read_trunc"));

        let subread_sets: HashMap<&str, HashSet<String>> = result
            .isoforms
            .iter()
            .map(|tx| (tx.name.as_str(), decoded_subreads(tx)))
            .collect();

        assert!(subread_sets.get("ref_a").unwrap().contains("read_trunc"));
        assert!(subread_sets.get("ref_b").unwrap().contains("read_trunc"));
    }

    #[test]
    fn sl_supported_same_5prime_single_exon_reads_still_merge() {
        let mut tracks = vec![
            make_track("long_sl", Strand::Plus, &[(100, 200)], "nanopore_read", 12),
            make_track("short_sl", Strand::Plus, &[(105, 180)], "nanopore_read", 12),
            make_track("same_sl", Strand::Plus, &[(100, 200)], "nanopore_read", 12),
        ];

        let keep = junction_simple_merge(&mut tracks, 11);

        assert_eq!(keep, vec![2]);
        assert!(has_subread(&tracks[2], "long_sl"));
        assert!(has_subread(&tracks[2], "short_sl"));
        assert!(has_subread(&tracks[2], "same_sl"));
    }

    #[test]
    fn sl_supported_same_junction_reads_merge_with_small_5prime_offset() {
        let mut tracks = vec![
            make_track(
                "long_sl",
                Strand::Plus,
                &[(100, 110), (120, 130), (140, 160)],
                "nanopore_read",
                12,
            ),
            make_track(
                "short_sl",
                Strand::Plus,
                &[(105, 110), (120, 130), (140, 150)],
                "nanopore_read",
                12,
            ),
        ];

        let keep = junction_simple_merge(&mut tracks, 11);

        assert_eq!(keep, vec![0]);
        assert!(has_subread(&tracks[0], "long_sl"));
        assert!(has_subread(&tracks[0], "short_sl"));
    }

    #[test]
    fn sl_supported_same_junction_singleton_large_5prime_offset_merges() {
        let mut tracks = vec![
            make_track(
                "long_sl",
                Strand::Plus,
                &[(50, 110), (120, 130), (140, 160)],
                "nanopore_read",
                12,
            ),
            make_track(
                "short_sl",
                Strand::Plus,
                &[(85, 110), (120, 130), (140, 150)],
                "nanopore_read",
                12,
            ),
        ];

        let keep = junction_simple_merge(&mut tracks, 11);

        assert_eq!(keep, vec![0]);
        assert!(has_subread(&tracks[0], "short_sl"));
    }

    #[test]
    fn sl_supported_same_junction_cluster_with_large_5prime_offset_is_retained() {
        let refs = vec![make_tx(
            "ref",
            Strand::Plus,
            &[(100, 110), (120, 130), (140, 150)],
            "isoform_anno",
            100,
        )];
        let reads = vec![
            make_tx(
                "alt_sl_a",
                Strand::Plus,
                &[(50, 110), (120, 130), (140, 150)],
                "nanopore_read",
                12,
            ),
            make_tx(
                "alt_sl_b",
                Strand::Plus,
                &[(52, 110), (120, 130), (140, 150)],
                "nanopore_read",
                12,
            ),
        ];

        let result = clusterj_with_name2_mode(&reads, Some(&refs), 1, 11, 0, 1, Name2Mode::Full);
        let iso_names: HashSet<_> = result.isoforms.iter().map(|tx| tx.name.as_str()).collect();

        assert!(iso_names.contains("ref"));
        let retained = mapped_isoform_id(&result, "alt_sl_a");
        assert_eq!(retained, mapped_isoform_id(&result, "alt_sl_b"));
        assert!(retained.starts_with(crate::identity::NOVEL_ISOFORM_PREFIX));
        assert!(iso_names.contains(retained));
    }

    #[test]
    fn supported_same_junction_three_prime_cluster_is_retained() {
        let refs = vec![make_tx(
            "ref",
            Strand::Minus,
            &[(100, 200), (300, 350), (400, 500)],
            "isoform_anno",
            100,
        )];
        let reads = vec![
            make_tx(
                "early_stop_a",
                Strand::Minus,
                &[(160, 200), (300, 350), (400, 500)],
                "nanopore_read",
                0,
            ),
            make_tx(
                "early_stop_b",
                Strand::Minus,
                &[(162, 200), (300, 350), (400, 500)],
                "nanopore_read",
                0,
            ),
            make_tx(
                "early_stop_c",
                Strand::Minus,
                &[(164, 200), (300, 350), (400, 500)],
                "nanopore_read",
                0,
            ),
            make_tx(
                "early_stop_d",
                Strand::Minus,
                &[(166, 200), (300, 350), (400, 500)],
                "nanopore_read",
                0,
            ),
            make_tx(
                "early_stop_e",
                Strand::Minus,
                &[(168, 200), (300, 350), (400, 500)],
                "nanopore_read",
                0,
            ),
        ];

        let result = clusterj_with_name2_mode(&reads, Some(&refs), 1, 11, 0, 1, Name2Mode::Full);
        let iso_names: HashSet<_> = result.isoforms.iter().map(|tx| tx.name.as_str()).collect();

        assert!(iso_names.contains("ref"));
        let retained = mapped_isoform_id(&result, "early_stop_a");
        assert!(retained.starts_with(crate::identity::NOVEL_ISOFORM_PREFIX));
        for read in [
            "early_stop_b",
            "early_stop_c",
            "early_stop_d",
            "early_stop_e",
        ] {
            assert_eq!(mapped_isoform_id(&result, read), retained);
        }

        let early_stop_tx = result
            .isoforms
            .iter()
            .find(|tx| tx.name == retained)
            .expect("early-stop isoform missing");
        let early_stop_reads = decoded_subreads(early_stop_tx);
        assert!(early_stop_reads.contains("early_stop_b"));
        assert!(early_stop_reads.contains("early_stop_e"));
    }

    #[test]
    fn same_junction_read_merges_into_reference_when_read_is_longer() {
        let refs = vec![make_tx(
            "ref",
            Strand::Plus,
            &[(100, 110), (120, 130), (140, 150)],
            "isoform_anno",
            100,
        )];
        let reads = vec![make_tx(
            "read_longer",
            Strand::Plus,
            &[(97, 110), (120, 130), (140, 153)],
            "nanopore_read",
            12,
        )];

        let result = clusterj_with_name2_mode(&reads, Some(&refs), 1, 11, 0, 1, Name2Mode::Full);
        let iso_names: HashSet<_> = result.isoforms.iter().map(|tx| tx.name.as_str()).collect();

        assert!(iso_names.contains("ref"));
        assert!(!iso_names.contains("read_longer"));

        let ref_tx = result
            .isoforms
            .iter()
            .find(|tx| tx.name == "ref")
            .expect("ref isoform missing");
        assert!(ref_tx
            .extra_fields
            .first()
            .is_some_and(|name2| name2.contains("read_longer")));
    }

    #[test]
    fn sl_supported_same_junction_read_merges_within_wide_5prime_offset() {
        let refs = vec![make_tx(
            "ref",
            Strand::Plus,
            &[(100, 110), (120, 130), (140, 150)],
            "isoform_anno",
            100,
        )];
        let reads = vec![make_tx(
            "read_19bp_offset",
            Strand::Plus,
            &[(81, 110), (120, 130), (140, 150)],
            "nanopore_read",
            12,
        )];

        let result = clusterj_with_name2_mode(&reads, Some(&refs), 1, 11, 0, 1, Name2Mode::Full);
        let iso_names: HashSet<_> = result.isoforms.iter().map(|tx| tx.name.as_str()).collect();

        assert!(iso_names.contains("ref"));
        assert!(!iso_names.contains("read_19bp_offset"));
    }

    #[test]
    fn sl_supported_singleton_5prime_junction_variant_merges_as_possible_degradation() {
        let refs = vec![make_tx(
            "ref",
            Strand::Plus,
            &[(100, 110), (120, 130), (140, 150)],
            "isoform_anno",
            100,
        )];
        let reads = vec![make_tx(
            "read_sl_trunc",
            Strand::Plus,
            &[(120, 130), (140, 150)],
            "nanopore_read",
            12,
        )];

        let result = clusterj_with_name2_mode(&reads, Some(&refs), 1, 11, 0, 1, Name2Mode::Full);
        let iso_names: HashSet<_> = result.isoforms.iter().map(|tx| tx.name.as_str()).collect();

        assert!(iso_names.contains("ref"));
        assert!(!iso_names.contains("read_sl_trunc"));
    }

    #[test]
    fn sl_supported_5prime_junction_variant_cluster_is_retained() {
        let refs = vec![make_tx(
            "ref",
            Strand::Plus,
            &[(100, 110), (120, 130), (140, 150)],
            "isoform_anno",
            100,
        )];
        let reads = vec![
            make_tx(
                "read_sl_trunc_a",
                Strand::Plus,
                &[(120, 130), (140, 150)],
                "nanopore_read",
                12,
            ),
            make_tx(
                "read_sl_trunc_b",
                Strand::Plus,
                &[(120, 130), (140, 150)],
                "nanopore_read",
                12,
            ),
        ];

        let result = clusterj_with_name2_mode(&reads, Some(&refs), 1, 11, 0, 1, Name2Mode::Full);
        let iso_names: HashSet<_> = result.isoforms.iter().map(|tx| tx.name.as_str()).collect();

        assert!(iso_names.contains("ref"));
        let retained = mapped_isoform_id(&result, "read_sl_trunc_a");
        assert_eq!(retained, mapped_isoform_id(&result, "read_sl_trunc_b"));
        assert!(retained.starts_with(crate::identity::NOVEL_ISOFORM_PREFIX));
    }

    #[test]
    fn batching_preserves_distinct_alignments_with_the_same_read_id() {
        let refs = vec![make_tx(
            "ref",
            Strand::Plus,
            &[(80, 170)],
            "isoform_anno",
            100,
        )];
        let reads = vec![
            make_tx(
                "duplicate_id",
                Strand::Plus,
                &[(50, 150)],
                "nanopore_read",
                100,
            ),
            make_tx(
                "duplicate_id",
                Strand::Plus,
                &[(100, 200)],
                "nanopore_read",
                100,
            ),
        ];

        let unbatched =
            clusterj_with_name2_mode(&reads, Some(&refs), 1, 11, 0, 10, Name2Mode::Full);
        let batched = clusterj_with_name2_mode(&reads, Some(&refs), 1, 11, 1, 10, Name2Mode::Full);

        assert_eq!(batched.isoforms, unbatched.isoforms);
        assert_eq!(batched.read_to_isoform, unbatched.read_to_isoform);
        let mapped = batched
            .read_to_isoform
            .iter()
            .filter(|(read_id, _)| read_id == "duplicate_id")
            .map(|(_, isoform_id)| isoform_id)
            .collect::<HashSet<_>>();
        assert_eq!(mapped.len(), 2);
    }

    #[test]
    fn sw_score_minus_one_removes_junction_sl_cluster_protection() {
        let refs = vec![make_tx(
            "ref",
            Strand::Plus,
            &[(100, 110), (120, 130), (140, 150)],
            "isoform_anno",
            100,
        )];
        let reads = vec![
            make_tx(
                "read_sl_trunc_a",
                Strand::Plus,
                &[(120, 130), (140, 150)],
                "nanopore_read",
                12,
            ),
            make_tx(
                "read_sl_trunc_b",
                Strand::Plus,
                &[(120, 130), (140, 150)],
                "nanopore_read",
                12,
            ),
        ];

        let protected = clusterj_with_name2_mode(&reads, Some(&refs), 1, 11, 0, 1, Name2Mode::Full);
        let protected_target = mapped_isoform_id(&protected, "read_sl_trunc_a");
        assert_eq!(
            protected_target,
            mapped_isoform_id(&protected, "read_sl_trunc_b")
        );
        assert!(protected_target.starts_with(crate::identity::NOVEL_ISOFORM_PREFIX));

        for batch_size in [0, 1] {
            let no_signal = clusterj_with_name2_mode(
                &reads,
                Some(&refs),
                1,
                -1,
                batch_size,
                10,
                Name2Mode::Full,
            );
            let no_signal_names: HashSet<_> = no_signal
                .isoforms
                .iter()
                .map(|tx| tx.name.as_str())
                .collect();
            assert!(no_signal_names.contains("ref"));
            assert!(!no_signal_names.contains("read_sl_trunc_a"));
            assert!(!no_signal_names.contains("read_sl_trunc_b"));
            assert!(no_signal
                .read_to_isoform
                .iter()
                .all(|(_, isoform_id)| isoform_id == "ref"));
        }
    }

    #[test]
    fn no_sl_default_merges_high_score_terminal_single_exon_reads() {
        let refs = vec![make_tx(
            "ref",
            Strand::Plus,
            &[(100, 150), (220, 260), (320, 350)],
            "isoform_anno",
            100,
        )];
        let reads = vec![
            make_tx(
                "terminal_a",
                Strand::Plus,
                &[(320, 350)],
                "nanopore_read",
                60,
            ),
            make_tx(
                "terminal_b",
                Strand::Plus,
                &[(332, 350)],
                "nanopore_read",
                60,
            ),
        ];

        let result = clusterj_with_name2_mode(
            &reads,
            Some(&refs),
            1,
            DEFAULT_SW_SCORE,
            0,
            1,
            Name2Mode::Full,
        );
        let iso_names: HashSet<_> = result.isoforms.iter().map(|tx| tx.name.as_str()).collect();

        assert!(iso_names.contains("ref"));
        assert!(!iso_names.contains("terminal_a"));
        assert!(!iso_names.contains("terminal_b"));
        assert_eq!(
            result.read_to_isoform,
            vec![
                ("terminal_a".to_owned(), "ref".to_owned()),
                ("terminal_b".to_owned(), "ref".to_owned()),
            ]
        );

        let counts =
            crate::count::count_by_read_to_isoform(&result.isoforms, &result.read_to_isoform)
                .unwrap();
        let ref_count = counts
            .iter()
            .find(|record| record.isoform_id == "ref")
            .expect("ref count missing");
        assert_eq!(ref_count.count, 2.0);
    }

    #[test]
    fn sl_supported_terminal_single_exon_cluster_is_retained() {
        let refs = vec![make_tx(
            "ref",
            Strand::Plus,
            &[(100, 150), (220, 260), (320, 350)],
            "isoform_anno",
            100,
        )];
        let reads = vec![
            make_tx(
                "terminal_sl_a",
                Strand::Plus,
                &[(320, 350)],
                "nanopore_read",
                12,
            ),
            make_tx(
                "terminal_sl_b",
                Strand::Plus,
                &[(320, 350)],
                "nanopore_read",
                12,
            ),
        ];

        let result = clusterj_with_name2_mode(&reads, Some(&refs), 1, 11, 0, 1, Name2Mode::Full);
        let iso_names: HashSet<_> = result.isoforms.iter().map(|tx| tx.name.as_str()).collect();

        assert!(iso_names.contains("ref"));
        let retained = mapped_isoform_id(&result, "terminal_sl_a");
        assert_eq!(retained, mapped_isoform_id(&result, "terminal_sl_b"));
        assert!(retained.starts_with(crate::identity::NOVEL_ISOFORM_PREFIX));
        assert_eq!(
            result.read_to_isoform,
            vec![
                ("terminal_sl_a".to_owned(), retained.to_owned()),
                ("terminal_sl_b".to_owned(), retained.to_owned()),
            ]
        );
    }

    #[test]
    fn sw_score_minus_one_keeps_normal_junction_merge() {
        let refs = vec![make_tx(
            "ref",
            Strand::Plus,
            &[(100, 110), (120, 130), (140, 150)],
            "isoform_anno",
            100,
        )];
        let reads = vec![
            make_tx(
                "read_trunc",
                Strand::Plus,
                &[(120, 130), (140, 150)],
                "nanopore_read",
                1,
            ),
            make_tx(
                "read_full",
                Strand::Plus,
                &[(100, 110), (120, 130), (140, 150)],
                "nanopore_read",
                1,
            ),
        ];

        for batch_size in [0, 1] {
            let result = clusterj_with_name2_mode(
                &reads,
                Some(&refs),
                1,
                -1,
                batch_size,
                10,
                Name2Mode::Full,
            );
            let iso_names: HashSet<_> = result.isoforms.iter().map(|tx| tx.name.as_str()).collect();

            assert!(iso_names.contains("ref"));
            assert!(!iso_names.contains("read_trunc"));
            assert!(!iso_names.contains("read_full"));
            assert_eq!(
                result.read_to_isoform,
                vec![
                    ("read_full".to_owned(), "ref".to_owned()),
                    ("read_trunc".to_owned(), "ref".to_owned()),
                ]
            );
        }
    }

    #[test]
    fn explicit_provenance_protects_plain_reference_and_maps_annotated_read() {
        let refs = vec![make_plain_tx_on(
            "chr1",
            "ref_plain",
            Strand::Plus,
            &[(100, 110), (120, 130), (140, 150)],
            100,
        )];
        // Biological annotation must not override the fact that this record came from the reads
        // input.
        let reads = vec![make_tx(
            "read_with_reference_annotation",
            Strand::Plus,
            &[(100, 110), (120, 130), (140, 150)],
            "isoform_anno",
            1,
        )];

        let (result, summary) = cluster_with_summary(&reads, Some(&refs));

        assert_eq!(result.isoforms.len(), 1);
        assert_eq!(result.isoforms[0].name, "ref_plain");
        assert_eq!(
            result.read_to_isoform,
            vec![(
                "read_with_reference_annotation".to_owned(),
                "ref_plain".to_owned()
            )]
        );
        assert!(result.unused.is_empty());
        assert_eq!(summary.input_reads, 1);
        assert_eq!(summary.represented_reads, 1);
        assert_eq!(summary.mapping_rows, 1);
    }

    #[test]
    fn read_id_identical_to_reference_id_is_not_suppressed() {
        let refs = vec![make_plain_tx_on(
            "chr1",
            "shared_id",
            Strand::Plus,
            &[(100, 110), (120, 130)],
            100,
        )];
        let reads = vec![make_plain_tx_on(
            "chr1",
            "shared_id",
            Strand::Plus,
            &[(100, 110), (120, 130)],
            1,
        )];

        let (result, summary) = cluster_with_summary(&reads, Some(&refs));

        assert_eq!(
            result.read_to_isoform,
            vec![("shared_id".to_owned(), "shared_id".to_owned())]
        );
        assert!(result.unused.is_empty());
        assert_eq!(summary.represented_reads, 1);
        let counts =
            crate::count::count_by_read_to_isoform(&result.isoforms, &result.read_to_isoform)
                .unwrap();
        assert_eq!(counts.len(), 1);
        assert_eq!(counts[0].isoform_id, "shared_id");
        assert_eq!(counts[0].count, 1.0);
    }

    #[test]
    fn unmatched_chromosome_strand_unknown_strand_and_locus_are_unused() {
        let refs = vec![make_plain_tx_on(
            "chr1",
            "ref",
            Strand::Plus,
            &[(100, 110), (120, 130)],
            100,
        )];
        let reads = vec![
            make_plain_tx_on(
                "chr1",
                "read_match",
                Strand::Plus,
                &[(100, 110), (120, 130)],
                1,
            ),
            make_plain_tx_on(
                "chr2",
                "read_wrong_chrom",
                Strand::Plus,
                &[(100, 110), (120, 130)],
                1,
            ),
            make_plain_tx_on(
                "chr1",
                "read_wrong_strand",
                Strand::Minus,
                &[(100, 110), (120, 130)],
                1,
            ),
            make_plain_tx_on(
                "chr1",
                "read_unknown_strand",
                Strand::Unknown,
                &[(100, 110), (120, 130)],
                1,
            ),
            make_plain_tx_on(
                "chr1",
                "read_disjoint",
                Strand::Plus,
                &[(500, 510), (520, 530)],
                1,
            ),
        ];

        let (result, summary) = cluster_with_summary(&reads, Some(&refs));
        let unused_names: HashSet<&str> = result.unused.iter().map(|tx| tx.name.as_str()).collect();

        assert_eq!(
            result.read_to_isoform,
            vec![("read_match".to_owned(), "ref".to_owned())]
        );
        assert_eq!(
            unused_names,
            HashSet::from([
                "read_wrong_chrom",
                "read_wrong_strand",
                "read_unknown_strand",
                "read_disjoint",
            ])
        );
        assert_eq!(summary.input_reads, 5);
        assert_eq!(summary.represented_reads, 1);
        assert_eq!(summary.rare_reads, 0);
        assert_eq!(summary.unmatched_reads, 4);
        assert_eq!(summary.unused_reads, 4);
    }

    #[test]
    fn absent_and_empty_reference_catalogs_return_every_read_as_unused() {
        let reads = vec![
            make_plain_tx_on("chr1", "read_a", Strand::Plus, &[(100, 110), (120, 130)], 1),
            make_plain_tx_on("chr2", "read_b", Strand::Unknown, &[(200, 210)], 1),
        ];

        let empty_references: Vec<Transcript> = Vec::new();
        for references in [None, Some(empty_references.as_slice())] {
            let (result, summary) = cluster_with_summary(&reads, references);
            assert!(result.isoforms.is_empty());
            assert!(result.read_to_isoform.is_empty());
            assert_eq!(result.unused, reads);
            assert_eq!(summary.input_reads, 2);
            assert_eq!(summary.represented_reads, 0);
            assert_eq!(summary.unmatched_reads, 2);
            assert_eq!(summary.unused_reads, 2);
        }
    }

    #[test]
    fn duplicate_read_ids_are_conserved_as_distinct_instances() {
        let refs = vec![make_plain_tx_on(
            "chr1",
            "ref",
            Strand::Plus,
            &[(100, 110), (120, 130)],
            100,
        )];
        let reads = vec![
            make_plain_tx_on(
                "chr1",
                "duplicate_id",
                Strand::Plus,
                &[(100, 110), (120, 130)],
                1,
            ),
            make_plain_tx_on(
                "chr2",
                "duplicate_id",
                Strand::Plus,
                &[(100, 110), (120, 130)],
                1,
            ),
        ];

        let (result, summary) = cluster_with_summary(&reads, Some(&refs));

        assert_eq!(
            result.read_to_isoform,
            vec![("duplicate_id".to_owned(), "ref".to_owned())]
        );
        assert_eq!(result.unused.len(), 1);
        assert_eq!(result.unused[0].name, "duplicate_id");
        assert_eq!(summary.input_reads, 2);
        assert_eq!(summary.represented_reads, 1);
        assert_eq!(summary.unmatched_reads, 1);
        assert_eq!(summary.unused_reads, 1);
    }

    #[test]
    fn duplicate_mapped_read_ids_produce_one_mapping_row_per_instance() {
        let refs = vec![make_plain_tx_on(
            "chr1",
            "ref",
            Strand::Plus,
            &[(100, 110), (120, 130)],
            100,
        )];
        let duplicate = make_plain_tx_on(
            "chr1",
            "duplicate_id",
            Strand::Plus,
            &[(100, 110), (120, 130)],
            1,
        );
        let reads = vec![duplicate.clone(), duplicate];

        let (result, summary) = cluster_with_summary(&reads, Some(&refs));

        assert_eq!(
            result.read_to_isoform,
            vec![
                ("duplicate_id".to_owned(), "ref".to_owned()),
                ("duplicate_id".to_owned(), "ref".to_owned()),
            ]
        );
        assert!(result.unused.is_empty());
        assert_eq!(summary.input_reads, 2);
        assert_eq!(summary.represented_reads, 2);
        assert_eq!(summary.mapping_rows, 2);
    }

    #[test]
    fn summary_distinguishes_rare_reads_from_unmatched_reads() {
        let refs = vec![make_plain_tx_on(
            "chr1",
            "ref",
            Strand::Plus,
            &[(100, 110), (200, 210)],
            100,
        )];
        let reads = vec![make_plain_tx_on(
            "chr1",
            "rare_read",
            Strand::Plus,
            &[(100, 150), (201, 210)],
            1,
        )];

        let (result, summary) = cluster_with_summary(&reads, Some(&refs));

        assert!(result.read_to_isoform.is_empty());
        assert_eq!(result.unused, reads);
        assert_eq!(summary.input_reads, 1);
        assert_eq!(summary.represented_reads, 0);
        assert_eq!(summary.rare_reads, 1);
        assert_eq!(summary.unmatched_reads, 0);
        assert_eq!(summary.unused_reads, 1);
    }

    #[test]
    fn reference_is_never_discarded_when_support_cutoff_exceeds_reference_weight() {
        let refs = vec![make_plain_tx_on(
            "chr1",
            "ref",
            Strand::Plus,
            &[(100, 110), (200, 210)],
            100,
        )];
        let (result, summary) = clusterj_with_options_and_summary(
            &[],
            Some(&refs),
            1,
            DEFAULT_SW_SCORE,
            0,
            1,
            Name2Mode::Full,
            SlMergeOptions::default(),
            ThreePrimeMergeOptions::default(),
            JunctionCorrectionOptions {
                min_support: 100,
                offset: DEFAULT_JUNCTION_CORRECTION_OFFSET,
            },
        );

        assert_eq!(result.isoforms.len(), 1);
        assert_eq!(result.isoforms[0].name, "ref");
        assert!(result.read_to_isoform.is_empty());
        assert!(result.unused.is_empty());
        assert_eq!(summary, JunctionClusterSummary::default());
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

        let single = clusterj_with_name2_mode(&reads, Some(&refs), 1, 11, 0, 1, Name2Mode::Full);
        let threaded = clusterj_with_name2_mode(&reads, Some(&refs), 4, 11, 0, 1, Name2Mode::Full);
        assert_eq!(single.isoforms, threaded.isoforms);
        assert_eq!(single.read_to_isoform, threaded.read_to_isoform);
        assert_eq!(single.unused, threaded.unused);
    }

    #[test]
    fn stable_novel_ids_and_payloads_are_input_order_independent() {
        let refs = vec![make_tx(
            "ref",
            Strand::Plus,
            &[(100, 200)],
            "isoform_anno",
            100,
        )];
        let reads = vec![
            make_tx("read,z|%", Strand::Plus, &[(50, 150)], "nanopore_read", 1),
            make_tx("read,a", Strand::Plus, &[(50, 150)], "nanopore_read", 1),
        ];

        let forward = clusterj_with_name2_mode(&reads, Some(&refs), 1, 11, 0, 1, Name2Mode::Full);
        let mut reversed_reads = reads.clone();
        reversed_reads.reverse();
        let reversed =
            clusterj_with_name2_mode(&reversed_reads, Some(&refs), 4, 11, 0, 1, Name2Mode::Full);

        assert_eq!(forward.isoforms, reversed.isoforms);
        assert_eq!(forward.read_to_isoform, reversed.read_to_isoform);
        let novel_id = mapped_isoform_id(&forward, "read,z|%");
        assert_eq!(novel_id, mapped_isoform_id(&forward, "read,a"));
        assert!(novel_id.starts_with(crate::identity::NOVEL_ISOFORM_PREFIX));
        let novel = forward
            .isoforms
            .iter()
            .find(|tx| tx.name == novel_id)
            .unwrap();
        assert_eq!(
            decoded_subreads(novel),
            HashSet::from(["read,z|%".to_owned(), "read,a".to_owned()])
        );
        assert!(novel.extra_fields[0].contains("%2C"));
        assert!(novel.extra_fields[0].contains("%7C"));
        assert!(novel.extra_fields[0].contains("%25"));
    }

    #[test]
    fn clustering_rejects_duplicate_or_reserved_reference_ids() {
        let duplicate = vec![
            make_tx("same", Strand::Plus, &[(0, 10)], "isoform_anno", 100),
            make_tx("same", Strand::Plus, &[(20, 30)], "isoform_anno", 100),
        ];
        let error = try_clusterj(&[], Some(&duplicate), 1, DEFAULT_SW_SCORE, 0, 1).unwrap_err();
        assert!(error.to_string().contains("duplicate reference isoform id"));

        let reserved = vec![make_tx(
            &format!("{}claimed", crate::identity::NOVEL_ISOFORM_PREFIX),
            Strand::Plus,
            &[(0, 10)],
            "isoform_anno",
            100,
        )];
        let error = try_clusterj(&[], Some(&reserved), 1, DEFAULT_SW_SCORE, 0, 1).unwrap_err();
        assert!(error
            .to_string()
            .contains("reserved novel-isoform namespace"));
    }

    #[test]
    fn rejects_zero_threads_and_support_at_library_boundary() {
        assert!(try_clusterj(&[], None, 0, DEFAULT_SW_SCORE, 0, 1).is_err());

        assert!(try_clusterj_with_options(
            &[],
            None,
            1,
            DEFAULT_SW_SCORE,
            0,
            1,
            Name2Mode::Full,
            SlMergeOptions {
                min_five_prime_cluster_support: 0,
                ..SlMergeOptions::default()
            },
            ThreePrimeMergeOptions::default(),
            JunctionCorrectionOptions::default(),
        )
        .is_err());

        assert!(try_clusterj_with_options(
            &[],
            None,
            1,
            DEFAULT_SW_SCORE,
            0,
            1,
            Name2Mode::Full,
            SlMergeOptions::default(),
            ThreePrimeMergeOptions::default(),
            JunctionCorrectionOptions {
                min_support: 0,
                ..JunctionCorrectionOptions::default()
            },
        )
        .is_err());
    }

    #[test]
    fn single_exon_index_candidates_cover_naive_merges_and_stay_sparse() {
        let sl = SlMergeOptions::default();
        let mut tracks = vec![make_track(
            "long",
            Strand::Plus,
            &[(0, 5000)],
            "nanopore_read",
            1,
        )];
        for index in 0..40 {
            let start = 100 + index * 120;
            tracks.push(make_track(
                &format!("short_{index}"),
                Strand::Plus,
                &[(start, start + 40)],
                "nanopore_read",
                1,
            ));
        }
        let junctions: Vec<Vec<u32>> = tracks
            .iter()
            .map(|track| junction_positions(&track.tx))
            .collect();
        let eligible = vec![true; tracks.len()];
        let index = SingleExonTargetIndex::new(&tracks, &junctions, &eligible);
        let mut candidate_count = 0usize;
        let mut naive_merge_count = 0usize;
        for i in 0..tracks.len() {
            let indexed: HashSet<_> = index.candidates(i, &tracks[i].tx, sl).into_iter().collect();
            candidate_count += indexed.len();
            for j in 0..tracks.len() {
                if i == j {
                    continue;
                }
                if single_exon_merge_kind(&tracks[i].tx, &tracks[j].tx, &junctions[j], sl).is_some()
                {
                    naive_merge_count += 1;
                    assert!(
                        indexed.contains(&j),
                        "index missed merge candidate {i}->{j}"
                    );
                }
            }
        }
        let naive_pairs = tracks.len() * (tracks.len() - 1);
        assert!(naive_merge_count > 0);
        assert!(
            candidate_count < naive_pairs / 4,
            "candidate_count={candidate_count} naive_pairs={naive_pairs}"
        );
    }

    #[test]
    fn fuzzy_same_junction_window_skips_far_first_junctions() {
        let junctions = [
            vec![100u32, 200, 300],
            vec![102, 204, 298],
            vec![400, 500, 600],
        ];
        let eligible = [true, true, true];
        let index = build_junction_length_index(&junctions, &eligible);
        let near = same_length_window_candidates(&index, &junctions[0], 5);
        assert!(near.contains(&1));
        assert!(!near.contains(&2));
        assert!(!junctions_equal(&junctions[0], &junctions[2], 5));
        assert!(junctions_equal(&junctions[0], &junctions[1], 5));
    }

    #[test]
    fn batched_merge_collapses_cross_batch_same_junction_reads() {
        let mut tracks = vec![make_track(
            "ref",
            Strand::Plus,
            &[(50, 400)],
            "isoform_anno",
            100,
        )];
        for index in 0..6 {
            tracks.push(make_track(
                &format!("read_{index}"),
                Strand::Plus,
                &[(100, 150), (200, 250), (300, 350)],
                "nanopore_read",
                1,
            ));
        }
        let merged = batch_junction_simple_merge(
            tracks,
            DEFAULT_SW_SCORE,
            2,
            100,
            SlMergeOptions::default(),
            ThreePrimeMergeOptions::default(),
            10,
        );
        let isoform = merged
            .iter()
            .find(|track| track.is_read())
            .expect("merged read isoform");
        assert_eq!(merged.iter().filter(|track| track.is_read()).count(), 1);
        for index in 0..6 {
            assert!(has_subread(isoform, &format!("read_{index}")));
        }
    }

    #[test]
    fn locus_read_cap_sends_dropped_reads_to_unused() {
        let refs = vec![make_tx(
            "ref",
            Strand::Plus,
            &[(0, 500)],
            "isoform_anno",
            100,
        )];
        let reads = (0..10)
            .map(|index| {
                make_tx(
                    &format!("read_{index}"),
                    Strand::Plus,
                    &[(100, 200)],
                    "nanopore_read",
                    1,
                )
            })
            .collect::<Vec<_>>();
        let (result, summary) = try_clusterj_with_runtime_options_and_summary(
            &reads,
            Some(&refs),
            1,
            DEFAULT_SW_SCORE,
            500,
            100,
            Name2Mode::Coverage,
            SlMergeOptions::default(),
            ThreePrimeMergeOptions::default(),
            JunctionCorrectionOptions::default(),
            ClusterjRuntimeOptions {
                max_reads_per_locus: 3,
                downsample_seed: 7,
                heartbeat_seconds: 0,
                heartbeat_top: 5,
            },
        )
        .unwrap();
        assert_eq!(summary.represented_reads, 3);
        assert_eq!(summary.unmatched_reads, 0);
        assert_eq!(summary.unused_reads, 7);
        assert_eq!(result.unused.len(), 7);
        assert_eq!(result.read_to_isoform.len(), 3);
    }
}
