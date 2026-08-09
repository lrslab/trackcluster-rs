//! Unique read-assignment aggregation of normalized modification observations.

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::io::Write;

use anyhow::Context;
use csv::WriterBuilder;

use crate::io::manifest::SampleRow;
use crate::io::mod_calls::ObservationReadResult;
use crate::model::{Coord, Strand, Transcript};
use crate::modification::{
    AssayMetadata, CoverageBasis, EligibilityReason, ImplicitSkipPolicy, ModSiteKey,
    ObservationState, SiteState,
};

/// Loaded normalized data for one `(sample, assay)` manifest row.
#[derive(Clone, Debug)]
pub struct ModSampleInput {
    /// Sample identifier.
    pub sample: String,
    /// Assay identifier.
    pub assay_id: String,
    /// Validated assay metadata.
    pub metadata: AssayMetadata,
    /// Parsed and exactly deduplicated observations.
    pub observations: ObservationReadResult,
    /// Optional exact primary-alignment coverage for this sample.
    pub coverage: Option<crate::io::coverage::BamCoverageResult>,
}

/// Scientific and QC gates for V1 modification aggregation.
#[derive(Clone, Debug)]
pub struct AggregateOptions {
    /// Per-assay hard-call probability thresholds.
    pub analysis_thresholds: BTreeMap<String, f64>,
    /// Minimum callable molecules required for an eligible site row.
    pub min_callable: u64,
    /// Minimum unique-read join rate for each sample/assay.
    pub min_read_join_rate: f64,
    /// Keep low-join rows as ineligible output instead of failing the command.
    pub allow_low_join: bool,
    /// Strand-oriented genomic bases keyed by every observed modification site.
    pub reference_bases: Option<BTreeMap<ModSiteKey, u8>>,
}

impl Default for AggregateOptions {
    fn default() -> Self {
        Self {
            analysis_thresholds: BTreeMap::new(),
            min_callable: 1,
            min_read_join_rate: 0.9,
            allow_low_join: false,
            reference_bases: None,
        }
    }
}

/// One sample/assay join audit row.
#[derive(Clone, Debug, PartialEq)]
pub struct JoinQcRow {
    /// Assay identifier.
    pub assay_id: String,
    /// Applied analysis threshold.
    pub analysis_threshold: f64,
    /// Sample identifier.
    pub sample: String,
    /// Input data rows, including exact duplicates.
    pub input_rows: usize,
    /// Unique normalized rows after exact deduplication.
    pub valid_rows: usize,
    /// Rows already represented in genomic coordinates.
    pub projected_rows: usize,
    /// Rows whose read ID joined a unique isoform assignment.
    pub joined_rows: usize,
    /// Distinct joined molecules.
    pub joined_reads: usize,
    /// Fraction of distinct input molecules that joined.
    pub read_join_rate: f64,
    /// Fraction of unique observation rows that joined.
    pub observation_join_rate: f64,
    /// Observation rows with no read assignment.
    pub unknown_read: usize,
    /// Rows rejected for an unknown sample (always zero after strict parsing).
    pub unknown_sample: usize,
    /// Rows whose mapping named a missing isoform (fail-fast, so zero on output).
    pub unknown_isoform: usize,
    /// Exact duplicate input rows folded by the normalized parser.
    pub duplicate_exact: usize,
    /// Conflicting duplicates (fail-fast, so zero on output).
    pub duplicate_conflict: usize,
    /// Source positions that could not be projected (zero for normalized genomic input).
    pub unprojectable: usize,
    /// Invalid probabilities (fail-fast, so zero on output).
    pub invalid_probability: usize,
    /// Whether the caller retained a complete candidate observation universe.
    pub candidate_observations_complete: bool,
}

/// Genomic-site observation-to-assignment join audit row.
#[derive(Clone, Debug, PartialEq)]
pub struct SiteJoinQcRow {
    /// Assay identifier.
    pub assay_id: String,
    /// Sample identifier.
    pub sample: String,
    /// Stable genomic site identifier.
    pub site_id: String,
    /// Genomic site/modification identity.
    pub site: ModSiteKey,
    /// Distinct normalized read-site observations before assignment join.
    pub input_rows: u64,
    /// Distinct normalized read-site observations with a unique isoform assignment.
    pub joined_rows: u64,
    /// Site-local observation join rate.
    pub observation_join_rate: f64,
    /// Whether the site-local join rate passes the configured threshold.
    pub passes_min_join_rate: bool,
}

/// Complete sample/isoform/site audit row.
#[derive(Clone, Debug, PartialEq)]
pub struct IsoformModSiteRow {
    /// Assay identifier.
    pub assay_id: String,
    /// Applied hard-call threshold.
    pub analysis_threshold: f64,
    /// Sample identifier.
    pub sample: String,
    /// Experimental group, if supplied.
    pub group: Option<String>,
    /// Gene identifier from the isoform catalog.
    pub gene: String,
    /// Isoform identifier.
    pub isoform_id: String,
    /// Stable genomic site identifier.
    pub site_id: String,
    /// Genomic site/modification identity.
    pub site: ModSiteKey,
    /// Candidate context or rule.
    pub context: Option<String>,
    /// Structural relationship between site and isoform.
    pub site_state: SiteState,
    /// Coverage evidence used for `n_covering`.
    pub coverage_basis: CoverageBasis,
    /// Molecules uniquely assigned to the sample/isoform.
    pub n_assigned: u64,
    /// Assigned molecules covering the genomic base, when available.
    pub n_covering: Option<u64>,
    /// Candidate observations evaluated by the caller.
    pub n_candidate: u64,
    /// Candidate observations with an interpretable hard-call state.
    pub n_callable: u64,
    /// Callable observations at or above the analysis threshold.
    pub n_modified: u64,
    /// Callable observations below the analysis threshold.
    pub n_unmodified: u64,
    /// Candidate observations excluded from the callable denominator.
    pub n_unknown: u64,
    /// Modified fraction when the candidate universe permits estimation.
    pub mod_fraction: Option<f64>,
    /// Mean explicit caller probability.
    pub mean_probability: Option<f64>,
    /// Wilson 95% interval lower bound.
    pub ci_low: Option<f64>,
    /// Wilson 95% interval upper bound.
    pub ci_high: Option<f64>,
    /// Stable eligibility reason; `ok` is the only eligible value.
    pub eligibility_reason: EligibilityReason,
}

/// Compact integer-count row for external statistical models.
#[derive(Clone, Debug, PartialEq)]
pub struct IsoformModDesignRow {
    /// Assay identifier.
    pub assay_id: String,
    /// Applied analysis threshold.
    pub analysis_threshold: f64,
    /// Sample identifier.
    pub sample: String,
    /// Experimental group, if supplied.
    pub group: Option<String>,
    /// Gene identifier.
    pub gene: String,
    /// Genomic site identifier.
    pub site_id: String,
    /// SAM-style modification code.
    pub mod_code: String,
    /// Isoform identifier.
    pub isoform_id: String,
    /// Modified molecule count.
    pub n_modified: u64,
    /// Unmodified molecule count.
    pub n_unmodified: u64,
    /// Modified fraction, or NA when ineligible for estimation.
    pub mod_fraction: Option<f64>,
    /// Stable eligibility reason.
    pub eligibility_reason: EligibilityReason,
}

/// All deterministic V1 aggregation outputs.
#[derive(Clone, Debug, Default)]
pub struct AggregateResult {
    /// Join audit rows.
    pub join_qc: Vec<JoinQcRow>,
    /// Site-local join audit rows, including sites with zero joined observations.
    pub site_join_qc: Vec<SiteJoinQcRow>,
    /// Complete sample/isoform/site rows.
    pub sites: Vec<IsoformModSiteRow>,
    /// Compact analysis-design rows.
    pub design: Vec<IsoformModDesignRow>,
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct SiteUniverseKey {
    assay_id: String,
    gene: String,
    site: ModSiteKey,
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct ObservationGroupKey {
    assay_id: String,
    sample: String,
    gene: String,
    isoform_id: String,
    site: ModSiteKey,
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct SiteJoinKey {
    assay_id: String,
    sample: String,
    site: ModSiteKey,
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct CoverageSiteIndexKey {
    assay_id: String,
    gene: String,
    chrom: String,
    strand: Strand,
}

#[derive(Clone, Debug, Default)]
struct ObservationCounts {
    explicit: u64,
    explicit_modified: u64,
    explicit_probability_sum: f64,
    implicit: u64,
    unknown: u64,
}

#[derive(Clone, Copy, Debug, Default)]
struct SiteJoinCounts {
    input_rows: u64,
    joined_rows: u64,
}

fn validate_options(options: &AggregateOptions) -> anyhow::Result<()> {
    if options.min_callable == 0 {
        anyhow::bail!("min_callable must be at least 1");
    }
    if !options.min_read_join_rate.is_finite() || !(0.0..=1.0).contains(&options.min_read_join_rate)
    {
        anyhow::bail!("min_read_join_rate must be finite and in [0, 1]");
    }
    for (assay_id, threshold) in &options.analysis_thresholds {
        if assay_id.trim().is_empty() || assay_id.chars().any(char::is_control) {
            anyhow::bail!("invalid analysis-threshold assay id {assay_id:?}");
        }
        if !threshold.is_finite() || !(0.0..=1.0).contains(threshold) {
            anyhow::bail!("analysis threshold for assay {assay_id:?} must be finite and in [0, 1]");
        }
    }
    Ok(())
}

fn isoform_gene(isoform: &Transcript) -> anyhow::Result<String> {
    let gene = isoform
        .metadata()
        .gene_id()
        .with_context(|| format!("isoform {:?} has no gene_id metadata", isoform.name))?;
    if gene.contains("||") {
        anyhow::bail!(
            "isoform {:?} has ambiguous multi-gene metadata {gene:?}",
            isoform.name
        );
    }
    Ok(gene.to_owned())
}

fn assay_compatible(left: &AssayMetadata, right: &AssayMetadata) -> bool {
    left.assay_id == right.assay_id
        && left.caller == right.caller
        && left.caller_version == right.caller_version
        && left.model_id == right.model_id
        && left.chemistry == right.chemistry
        && left.candidate_rule == right.candidate_rule
        && left.source_emission_threshold == right.source_emission_threshold
        && left.source_site_filter == right.source_site_filter
        && left.coordinate_source == right.coordinate_source
        && left.read_id_mapping == right.read_id_mapping
}

fn context_flank(context: Option<&str>, candidate_rule: &str) -> u32 {
    let context = context.unwrap_or(candidate_rule);
    if context.eq_ignore_ascii_case("DRACH") {
        return 2;
    }
    let bytes = context.as_bytes();
    if bytes.len() % 2 == 1
        && bytes.len() > 1
        && bytes.iter().all(|base| {
            matches!(
                base.to_ascii_uppercase(),
                b'A' | b'C'
                    | b'G'
                    | b'T'
                    | b'U'
                    | b'R'
                    | b'Y'
                    | b'S'
                    | b'W'
                    | b'K'
                    | b'M'
                    | b'B'
                    | b'D'
                    | b'H'
                    | b'V'
                    | b'N'
            )
        })
    {
        u32::try_from(bytes.len() / 2).unwrap_or(u32::MAX)
    } else {
        0
    }
}

fn classify_site(isoform: &Transcript, site: &ModSiteKey, flank: u32) -> SiteState {
    if isoform.strand == Strand::Unknown {
        return SiteState::Unprojectable;
    }
    if isoform.chrom != site.chrom || isoform.strand != site.strand {
        return SiteState::StructurallyAbsent;
    }
    let pos = Coord::new(site.pos0);
    let Some(exon) = isoform
        .exons
        .iter()
        .find(|exon| exon.start <= pos && pos < exon.end)
    else {
        return SiteState::StructurallyAbsent;
    };
    let left = site.pos0 - exon.start.get();
    let right = exon.end.get() - site.pos0 - 1;
    if left < flank || right < flank {
        SiteState::ContextDependent
    } else {
        SiteState::Present
    }
}

fn reference_base_matches_mod_code(reference_base: u8, mod_code: &str) -> bool {
    let canonical_base = mod_code
        .as_bytes()
        .first()
        .copied()
        .expect("validated modification code")
        .to_ascii_uppercase();
    canonical_base == b'N'
        || reference_base.to_ascii_uppercase()
            == if canonical_base == b'U' {
                b'T'
            } else {
                canonical_base
            }
}

fn wilson_interval(modified: u64, callable: u64) -> Option<(f64, f64)> {
    if callable == 0 {
        return None;
    }
    const Z: f64 = 1.959_963_984_540_054;
    let n = callable as f64;
    let fraction = modified as f64 / n;
    let z2 = Z * Z;
    let denominator = 1.0 + z2 / n;
    let center = (fraction + z2 / (2.0 * n)) / denominator;
    let margin = Z * ((fraction * (1.0 - fraction) / n + z2 / (4.0 * n * n)).sqrt()) / denominator;
    Some((
        (center - margin).max(0.0).min(fraction),
        (center + margin).min(1.0).max(fraction),
    ))
}

/// Aggregate normalized observations using one unique isoform assignment per molecule.
pub fn aggregate_modifications(
    samples: &[SampleRow],
    isoforms: &[Transcript],
    read_to_isoform: &[(String, String)],
    mod_inputs: &[ModSampleInput],
    options: &AggregateOptions,
) -> anyhow::Result<AggregateResult> {
    validate_options(options)?;
    if samples.is_empty() {
        anyhow::bail!("sample manifest is empty");
    }
    if isoforms.is_empty() {
        anyhow::bail!("isoform catalog is empty");
    }
    if mod_inputs.is_empty() {
        anyhow::bail!("modification manifest is empty");
    }

    let samples_by_id = samples
        .iter()
        .map(|sample| (sample.sample.as_str(), sample))
        .collect::<HashMap<_, _>>();
    if samples_by_id.len() != samples.len() {
        anyhow::bail!("sample manifest contains duplicate sample identifiers");
    }

    let mut isoforms_by_id: HashMap<&str, (&Transcript, String)> = HashMap::new();
    let mut isoforms_by_gene: BTreeMap<String, Vec<&Transcript>> = BTreeMap::new();
    for isoform in isoforms {
        let gene = isoform_gene(isoform)?;
        if isoforms_by_id
            .insert(isoform.name.as_str(), (isoform, gene.clone()))
            .is_some()
        {
            anyhow::bail!("duplicate isoform id {:?}", isoform.name);
        }
        isoforms_by_gene.entry(gene).or_default().push(isoform);
    }
    for values in isoforms_by_gene.values_mut() {
        values.sort_by(|left, right| left.name.cmp(&right.name));
    }

    let mut assignments: HashMap<String, String> = HashMap::new();
    let mut n_assigned: HashMap<(String, String), u64> = HashMap::new();
    for (read_id, isoform_id) in read_to_isoform {
        let (sample, _) = crate::sample::split_tagged_read_name(read_id).with_context(|| {
            format!("read-to-isoform id {read_id:?} must use <sample>::<read_id>")
        })?;
        if !samples_by_id.contains_key(sample) {
            anyhow::bail!("read-to-isoform id {read_id:?} has unknown sample prefix {sample:?}");
        }
        if !isoforms_by_id.contains_key(isoform_id.as_str()) {
            anyhow::bail!(
                "read-to-isoform id {read_id:?} references unknown isoform {isoform_id:?}"
            );
        }
        match assignments.get(read_id) {
            Some(existing) if existing != isoform_id => {
                anyhow::bail!(
                    "read {read_id:?} maps to multiple isoforms ({existing:?}, {isoform_id:?}); modification aggregation requires unique assignments"
                );
            }
            Some(_) => continue,
            None => {
                assignments.insert(read_id.clone(), isoform_id.clone());
                *n_assigned
                    .entry((sample.to_owned(), isoform_id.clone()))
                    .or_default() += 1;
            }
        }
    }

    let expected_assays = mod_inputs
        .iter()
        .map(|input| input.assay_id.clone())
        .collect::<BTreeSet<_>>();
    let threshold_assays = options
        .analysis_thresholds
        .keys()
        .cloned()
        .collect::<BTreeSet<_>>();
    if expected_assays != threshold_assays {
        anyhow::bail!(
            "analysis thresholds must match modification assays exactly; inputs={expected_assays:?}, thresholds={threshold_assays:?}"
        );
    }

    let mut seen_sample_assays = HashSet::new();
    let mut representative_metadata: HashMap<&str, &AssayMetadata> = HashMap::new();
    let mut metadata_by_sample_assay: HashMap<(String, String), &AssayMetadata> = HashMap::new();
    let mut assay_samples: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    for input in mod_inputs {
        input.metadata.validate().map_err(anyhow::Error::msg)?;
        if input.sample.trim().is_empty() || input.assay_id.trim().is_empty() {
            anyhow::bail!("modification sample and assay identifiers must not be empty");
        }
        if !samples_by_id.contains_key(input.sample.as_str()) {
            anyhow::bail!("modification input has unknown sample {:?}", input.sample);
        }
        if input.metadata.assay_id != input.assay_id {
            anyhow::bail!(
                "modification input assay {:?} does not match metadata assay {:?}",
                input.assay_id,
                input.metadata.assay_id
            );
        }
        let analysis_threshold = options.analysis_thresholds[&input.assay_id];
        if input
            .metadata
            .source_emission_threshold
            .is_some_and(|emission| analysis_threshold < emission)
        {
            anyhow::bail!(
                "analysis threshold {analysis_threshold} for assay {:?} is below source emission threshold {}; omitted candidates make hard-call counts incomplete",
                input.assay_id,
                input.metadata.source_emission_threshold.expect("checked above")
            );
        }
        if !seen_sample_assays.insert((input.sample.clone(), input.assay_id.clone())) {
            anyhow::bail!(
                "duplicate modification input for sample {:?}, assay {:?}",
                input.sample,
                input.assay_id
            );
        }
        if let Some(existing) = representative_metadata.get(input.assay_id.as_str()) {
            if !assay_compatible(existing, &input.metadata) {
                anyhow::bail!(
                    "incompatible metadata share assay_id {:?}; caller/model/chemistry/probability semantics must match",
                    input.assay_id
                );
            }
        } else {
            representative_metadata.insert(input.assay_id.as_str(), &input.metadata);
        }
        metadata_by_sample_assay.insert(
            (input.sample.clone(), input.assay_id.clone()),
            &input.metadata,
        );
        assay_samples
            .entry(input.assay_id.clone())
            .or_default()
            .insert(input.sample.clone());
    }

    let mut join_qc = Vec::new();
    let mut low_join = HashSet::new();
    let mut universe: BTreeMap<SiteUniverseKey, Option<String>> = BTreeMap::new();
    let mut observation_groups: BTreeMap<ObservationGroupKey, ObservationCounts> = BTreeMap::new();
    let mut site_join_counts: BTreeMap<SiteJoinKey, SiteJoinCounts> = BTreeMap::new();

    for input in mod_inputs {
        let threshold = options.analysis_thresholds[&input.assay_id];
        let mut input_reads = HashSet::new();
        let mut joined_reads = HashSet::new();
        let mut joined_rows = 0usize;
        let mut unknown_read = 0usize;

        for observation in &input.observations.observations {
            if observation.key.sample != input.sample || observation.key.assay_id != input.assay_id
            {
                anyhow::bail!(
                    "observation sample/assay ({:?}, {:?}) does not match mod manifest ({:?}, {:?})",
                    observation.key.sample,
                    observation.key.assay_id,
                    input.sample,
                    input.assay_id
                );
            }
            input_reads.insert(observation.key.read_id.as_str());
            let site_join = site_join_counts
                .entry(SiteJoinKey {
                    assay_id: input.assay_id.clone(),
                    sample: input.sample.clone(),
                    site: observation.key.site.clone(),
                })
                .or_default();
            site_join.input_rows += 1;
            let Some(isoform_id) = assignments.get(&observation.key.read_id) else {
                unknown_read += 1;
                continue;
            };
            site_join.joined_rows += 1;
            let (_, gene) = isoforms_by_id
                .get(isoform_id.as_str())
                .expect("assignment isoforms validated above");
            joined_rows += 1;
            joined_reads.insert(observation.key.read_id.as_str());

            let universe_key = SiteUniverseKey {
                assay_id: input.assay_id.clone(),
                gene: gene.clone(),
                site: observation.key.site.clone(),
            };
            match universe.get(&universe_key) {
                Some(existing) if existing != &observation.context => {
                    anyhow::bail!(
                        "conflicting context at assay/gene/site {:?}: {:?} vs {:?}",
                        universe_key,
                        existing,
                        observation.context
                    );
                }
                Some(_) => {}
                None => {
                    universe.insert(universe_key, observation.context.clone());
                }
            }

            let group = observation_groups
                .entry(ObservationGroupKey {
                    assay_id: input.assay_id.clone(),
                    sample: input.sample.clone(),
                    gene: gene.clone(),
                    isoform_id: isoform_id.clone(),
                    site: observation.key.site.clone(),
                })
                .or_default();
            match observation.observation_state {
                ObservationState::ExplicitProbability => {
                    let probability = observation
                        .probability
                        .expect("validated explicit probability");
                    group.explicit += 1;
                    group.explicit_probability_sum += probability;
                    if probability >= threshold {
                        group.explicit_modified += 1;
                    }
                }
                ObservationState::ImplicitBelowEmissionThreshold => group.implicit += 1,
                ObservationState::Unknown => group.unknown += 1,
            }
        }

        let read_join_rate = if input_reads.is_empty() {
            0.0
        } else {
            joined_reads.len() as f64 / input_reads.len() as f64
        };
        let observation_join_rate = if input.observations.observations.is_empty() {
            0.0
        } else {
            joined_rows as f64 / input.observations.observations.len() as f64
        };
        if read_join_rate < options.min_read_join_rate {
            low_join.insert((input.sample.clone(), input.assay_id.clone()));
            if !options.allow_low_join {
                anyhow::bail!(
                    "sample {:?}, assay {:?} read join rate {:.6} is below minimum {:.6}; use --allow-low-join to emit ineligible rows",
                    input.sample,
                    input.assay_id,
                    read_join_rate,
                    options.min_read_join_rate
                );
            }
        }
        join_qc.push(JoinQcRow {
            assay_id: input.assay_id.clone(),
            analysis_threshold: threshold,
            sample: input.sample.clone(),
            input_rows: input.observations.input_rows,
            valid_rows: input.observations.observations.len(),
            projected_rows: input.observations.observations.len(),
            joined_rows,
            joined_reads: joined_reads.len(),
            read_join_rate,
            observation_join_rate,
            unknown_read,
            unknown_sample: 0,
            unknown_isoform: 0,
            duplicate_exact: input.observations.duplicate_exact,
            duplicate_conflict: 0,
            unprojectable: 0,
            invalid_probability: 0,
            candidate_observations_complete: input.metadata.candidate_observations_complete,
        });
    }
    join_qc.sort_by(|left, right| {
        left.assay_id
            .cmp(&right.assay_id)
            .then_with(|| left.sample.cmp(&right.sample))
    });
    let site_join_qc = site_join_counts
        .iter()
        .map(|(key, counts)| {
            let observation_join_rate = counts.joined_rows as f64 / counts.input_rows as f64;
            SiteJoinQcRow {
                assay_id: key.assay_id.clone(),
                sample: key.sample.clone(),
                site_id: key.site.site_id(),
                site: key.site.clone(),
                input_rows: counts.input_rows,
                joined_rows: counts.joined_rows,
                observation_join_rate,
                passes_min_join_rate: observation_join_rate >= options.min_read_join_rate,
            }
        })
        .collect::<Vec<_>>();

    let mut exact_coverage_inputs = HashSet::new();
    let mut coverage_counts: HashMap<ObservationGroupKey, u64> = HashMap::new();
    let mut coverage_site_index: BTreeMap<CoverageSiteIndexKey, Vec<ModSiteKey>> = BTreeMap::new();
    for key in universe.keys() {
        coverage_site_index
            .entry(CoverageSiteIndexKey {
                assay_id: key.assay_id.clone(),
                gene: key.gene.clone(),
                chrom: key.site.chrom.clone(),
                strand: key.site.strand,
            })
            .or_default()
            .push(key.site.clone());
    }
    for sites in coverage_site_index.values_mut() {
        sites.sort();
    }

    for input in mod_inputs {
        let Some(coverage) = input.coverage.as_ref() else {
            continue;
        };
        exact_coverage_inputs.insert((input.sample.clone(), input.assay_id.clone()));
        let mut coverage_read_ids = HashSet::new();
        for read in &coverage.reads {
            let (sample, _) =
                crate::sample::split_tagged_read_name(&read.read_id).with_context(|| {
                    format!(
                        "coverage read id {:?} must use <sample>::<read_id>",
                        read.read_id
                    )
                })?;
            if sample != input.sample {
                anyhow::bail!(
                    "coverage read id {:?} has sample prefix {sample:?}, expected {:?}",
                    read.read_id,
                    input.sample
                );
            }
            if !coverage_read_ids.insert(read.read_id.as_str()) {
                anyhow::bail!("duplicate coverage read id {:?}", read.read_id);
            }
            if read.strand == Strand::Unknown || read.match_blocks.is_empty() {
                anyhow::bail!(
                    "coverage read {:?} must have a known strand and at least one match block",
                    read.read_id
                );
            }
            if read
                .match_blocks
                .windows(2)
                .any(|pair| pair[0].end > pair[1].start)
            {
                anyhow::bail!(
                    "coverage read {:?} has overlapping match blocks",
                    read.read_id
                );
            }
        }

        let mut missing_assignments = 0usize;
        let mut first_missing = None;
        for read_id in assignments.keys() {
            let (sample, _) = crate::sample::split_tagged_read_name(read_id)
                .expect("assignment read IDs validated above");
            if sample == input.sample && !coverage_read_ids.contains(read_id.as_str()) {
                missing_assignments += 1;
                first_missing.get_or_insert_with(|| read_id.clone());
            }
        }
        if missing_assignments > 0 {
            anyhow::bail!(
                "coverage BAM for sample {:?} is missing {missing_assignments} uniquely assigned reads; first missing read is {:?}",
                input.sample,
                first_missing.expect("count is positive")
            );
        }

        for read in &coverage.reads {
            let Some(isoform_id) = assignments.get(&read.read_id) else {
                continue;
            };
            let (isoform, gene) = isoforms_by_id
                .get(isoform_id.as_str())
                .expect("assignment isoforms validated above");
            if read.chrom != isoform.chrom || read.strand != isoform.strand {
                anyhow::bail!(
                    "coverage alignment for read {:?} is {}:{} but assigned isoform {:?} is {}:{}",
                    read.read_id,
                    read.chrom,
                    read.strand.as_char(),
                    isoform_id,
                    isoform.chrom,
                    isoform.strand.as_char()
                );
            }
            let index_key = CoverageSiteIndexKey {
                assay_id: input.assay_id.clone(),
                gene: gene.clone(),
                chrom: read.chrom.clone(),
                strand: read.strand,
            };
            let Some(sites) = coverage_site_index.get(&index_key) else {
                continue;
            };
            for block in &read.match_blocks {
                let start = sites.partition_point(|site| site.pos0 < block.start.get());
                let end = sites.partition_point(|site| site.pos0 < block.end.get());
                for site in &sites[start..end] {
                    *coverage_counts
                        .entry(ObservationGroupKey {
                            assay_id: input.assay_id.clone(),
                            sample: input.sample.clone(),
                            gene: gene.clone(),
                            isoform_id: isoform_id.clone(),
                            site: site.clone(),
                        })
                        .or_default() += 1;
                }
            }
        }
    }

    let mut sites = Vec::new();
    for (universe_key, context) in universe {
        let metadata = representative_metadata[universe_key.assay_id.as_str()];
        let threshold = options.analysis_thresholds[&universe_key.assay_id];
        let flank = context_flank(context.as_deref(), &metadata.candidate_rule);
        let gene_isoforms = &isoforms_by_gene[&universe_key.gene];
        for sample in &assay_samples[&universe_key.assay_id] {
            let sample_row = samples_by_id[sample.as_str()];
            let sample_metadata =
                metadata_by_sample_assay[&(sample.clone(), universe_key.assay_id.clone())];
            for isoform in gene_isoforms {
                let group_key = ObservationGroupKey {
                    assay_id: universe_key.assay_id.clone(),
                    sample: sample.clone(),
                    gene: universe_key.gene.clone(),
                    isoform_id: isoform.name.clone(),
                    site: universe_key.site.clone(),
                };
                let counts = observation_groups.get(&group_key);
                let structural_state = classify_site(isoform, &universe_key.site, flank);
                let reference_mismatch = options
                    .reference_bases
                    .as_ref()
                    .map(|bases| {
                        bases
                            .get(&universe_key.site)
                            .copied()
                            .with_context(|| {
                                format!(
                                    "reference base was not loaded for observed site {:?}",
                                    universe_key.site
                                )
                            })
                            .map(|base| {
                                !reference_base_matches_mod_code(base, &universe_key.site.mod_code)
                            })
                    })
                    .transpose()?
                    .unwrap_or(false);
                let site_state = if matches!(
                    structural_state,
                    SiteState::Present | SiteState::ContextDependent
                ) && reference_mismatch
                {
                    SiteState::ReferenceBaseMismatch
                } else {
                    structural_state
                };
                let explicit = counts.map_or(0, |value| value.explicit);
                let explicit_modified = counts.map_or(0, |value| value.explicit_modified);
                let explicit_probability_sum =
                    counts.map_or(0.0, |value| value.explicit_probability_sum);
                let implicit = counts.map_or(0, |value| value.implicit);
                let source_unknown = counts.map_or(0, |value| value.unknown);
                let n_candidate = explicit + implicit + source_unknown;
                let has_exact_coverage = exact_coverage_inputs
                    .contains(&(sample.clone(), universe_key.assay_id.clone()));
                let n_covering = has_exact_coverage
                    .then(|| coverage_counts.get(&group_key).copied().unwrap_or_default());
                if n_covering.is_some_and(|covering| n_candidate > covering) {
                    anyhow::bail!(
                        "candidate observations exceed exact coverage for sample {:?}, isoform {:?}, site {:?}: candidates={n_candidate}, covering={}",
                        sample,
                        isoform.name,
                        universe_key.site,
                        n_covering.expect("checked above")
                    );
                }

                let site_can_be_called =
                    matches!(site_state, SiteState::Present | SiteState::ContextDependent);
                let implicit_callable = site_can_be_called
                    && sample_metadata.implicit_skip_policy == ImplicitSkipPolicy::LowProbability
                    && sample_metadata
                        .source_emission_threshold
                        .is_some_and(|emission| threshold >= emission);
                let (n_modified, explicit_unmodified) = if site_can_be_called {
                    (explicit_modified, explicit - explicit_modified)
                } else {
                    (0, 0)
                };
                let n_unmodified =
                    explicit_unmodified + if implicit_callable { implicit } else { 0 };
                let n_callable = n_modified + n_unmodified;
                let n_unknown = n_candidate - n_callable;
                debug_assert_eq!(n_candidate, n_callable + n_unknown);

                let mean_probability = if site_can_be_called && explicit > 0 {
                    Some(explicit_probability_sum / explicit as f64)
                } else {
                    None
                };
                let join_is_low =
                    low_join.contains(&(sample.clone(), universe_key.assay_id.clone()));
                let site_join_is_low = site_join_counts
                    .get(&SiteJoinKey {
                        assay_id: universe_key.assay_id.clone(),
                        sample: sample.clone(),
                        site: universe_key.site.clone(),
                    })
                    .is_some_and(|counts| {
                        counts.joined_rows as f64 / (counts.input_rows as f64)
                            < options.min_read_join_rate
                    });
                let denominator_is_known = source_unknown == 0
                    && (implicit == 0 || implicit_callable)
                    && sample_metadata.implicit_skip_policy != ImplicitSkipPolicy::Unknown;
                let reason = match site_state {
                    SiteState::StructurallyAbsent => EligibilityReason::SiteAbsent,
                    SiteState::ContextDependent => EligibilityReason::ContextDependent,
                    SiteState::ReferenceBaseMismatch => EligibilityReason::ReferenceBaseMismatch,
                    SiteState::Unprojectable => EligibilityReason::Unprojectable,
                    SiteState::Present if !sample_metadata.candidate_observations_complete => {
                        EligibilityReason::IncompleteCandidateUniverse
                    }
                    SiteState::Present if join_is_low => EligibilityReason::JoinRateLow,
                    SiteState::Present if site_join_is_low => EligibilityReason::SiteJoinRateLow,
                    SiteState::Present if !denominator_is_known => {
                        EligibilityReason::UnknownDenominator
                    }
                    SiteState::Present if n_callable < options.min_callable => {
                        EligibilityReason::LowCallable
                    }
                    SiteState::Present => EligibilityReason::Ok,
                };
                let fraction_is_defined = n_callable > 0
                    && sample_metadata.candidate_observations_complete
                    && denominator_is_known
                    && site_can_be_called;
                let mod_fraction =
                    fraction_is_defined.then(|| n_modified as f64 / n_callable as f64);
                let (ci_low, ci_high) = if fraction_is_defined {
                    wilson_interval(n_modified, n_callable)
                        .map(|(low, high)| (Some(low), Some(high)))
                        .unwrap_or((None, None))
                } else {
                    (None, None)
                };
                sites.push(IsoformModSiteRow {
                    assay_id: universe_key.assay_id.clone(),
                    analysis_threshold: threshold,
                    sample: sample.clone(),
                    group: sample_row.group.clone(),
                    gene: universe_key.gene.clone(),
                    isoform_id: isoform.name.clone(),
                    site_id: universe_key.site.site_id(),
                    site: universe_key.site.clone(),
                    context: context.clone(),
                    site_state,
                    coverage_basis: if has_exact_coverage {
                        CoverageBasis::BamExact
                    } else {
                        CoverageBasis::Unavailable
                    },
                    n_assigned: n_assigned
                        .get(&(sample.clone(), isoform.name.clone()))
                        .copied()
                        .unwrap_or(0),
                    n_covering,
                    n_candidate,
                    n_callable,
                    n_modified,
                    n_unmodified,
                    n_unknown,
                    mod_fraction,
                    mean_probability,
                    ci_low,
                    ci_high,
                    eligibility_reason: reason,
                });
            }
        }
    }
    sites.sort_by(|left, right| {
        left.assay_id
            .cmp(&right.assay_id)
            .then_with(|| left.sample.cmp(&right.sample))
            .then_with(|| left.gene.cmp(&right.gene))
            .then_with(|| left.isoform_id.cmp(&right.isoform_id))
            .then_with(|| left.site.cmp(&right.site))
    });
    let design = sites
        .iter()
        .map(|row| IsoformModDesignRow {
            assay_id: row.assay_id.clone(),
            analysis_threshold: row.analysis_threshold,
            sample: row.sample.clone(),
            group: row.group.clone(),
            gene: row.gene.clone(),
            site_id: row.site_id.clone(),
            mod_code: row.site.mod_code.clone(),
            isoform_id: row.isoform_id.clone(),
            n_modified: row.n_modified,
            n_unmodified: row.n_unmodified,
            mod_fraction: row.mod_fraction,
            eligibility_reason: row.eligibility_reason,
        })
        .collect();

    Ok(AggregateResult {
        join_qc,
        site_join_qc,
        sites,
        design,
    })
}

fn optional<T: ToString>(value: Option<T>) -> String {
    value
        .map(|value| value.to_string())
        .unwrap_or_else(|| "NA".to_owned())
}

/// Write join/QC rows as deterministic TSV.
pub fn write_join_qc_tsv<W: Write>(writer: W, rows: &[JoinQcRow]) -> anyhow::Result<()> {
    let mut output = WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_writer(writer);
    output.write_record([
        "assay_id",
        "analysis_threshold",
        "sample",
        "input_rows",
        "valid_rows",
        "projected_rows",
        "joined_rows",
        "joined_reads",
        "read_join_rate",
        "observation_join_rate",
        "unknown_read",
        "unknown_sample",
        "unknown_isoform",
        "duplicate_exact",
        "duplicate_conflict",
        "unprojectable",
        "invalid_probability",
        "candidate_observations_complete",
    ])?;
    for row in rows {
        output.write_record([
            row.assay_id.clone(),
            row.analysis_threshold.to_string(),
            row.sample.clone(),
            row.input_rows.to_string(),
            row.valid_rows.to_string(),
            row.projected_rows.to_string(),
            row.joined_rows.to_string(),
            row.joined_reads.to_string(),
            row.read_join_rate.to_string(),
            row.observation_join_rate.to_string(),
            row.unknown_read.to_string(),
            row.unknown_sample.to_string(),
            row.unknown_isoform.to_string(),
            row.duplicate_exact.to_string(),
            row.duplicate_conflict.to_string(),
            row.unprojectable.to_string(),
            row.invalid_probability.to_string(),
            row.candidate_observations_complete.to_string(),
        ])?;
    }
    output.flush().context("flush modification join QC")?;
    Ok(())
}

/// Write site-local join/QC rows as deterministic TSV.
pub fn write_site_join_qc_tsv<W: Write>(writer: W, rows: &[SiteJoinQcRow]) -> anyhow::Result<()> {
    let mut output = WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_writer(writer);
    output.write_record([
        "assay_id",
        "sample",
        "site_id",
        "chrom",
        "pos0",
        "strand",
        "mod_code",
        "input_rows",
        "joined_rows",
        "observation_join_rate",
        "passes_min_join_rate",
    ])?;
    for row in rows {
        output.write_record([
            row.assay_id.clone(),
            row.sample.clone(),
            row.site_id.clone(),
            row.site.chrom.clone(),
            row.site.pos0.to_string(),
            row.site.strand.as_char().to_string(),
            row.site.mod_code.clone(),
            row.input_rows.to_string(),
            row.joined_rows.to_string(),
            row.observation_join_rate.to_string(),
            row.passes_min_join_rate.to_string(),
        ])?;
    }
    output
        .flush()
        .context("flush modification site-local join QC")?;
    Ok(())
}

/// Write complete sample/isoform/site rows as deterministic TSV.
pub fn write_isoform_mod_sites_tsv<W: Write>(
    writer: W,
    rows: &[IsoformModSiteRow],
) -> anyhow::Result<()> {
    let mut output = WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_writer(writer);
    output.write_record([
        "assay_id",
        "analysis_threshold",
        "sample",
        "group",
        "gene",
        "isoform_id",
        "site_id",
        "chrom",
        "pos0",
        "strand",
        "mod_code",
        "context",
        "site_state",
        "coverage_basis",
        "n_assigned",
        "n_covering",
        "n_candidate",
        "n_callable",
        "n_modified",
        "n_unmodified",
        "n_unknown",
        "mod_fraction",
        "mean_probability",
        "ci_low",
        "ci_high",
        "eligibility",
        "eligibility_reason",
    ])?;
    for row in rows {
        output.write_record([
            row.assay_id.clone(),
            row.analysis_threshold.to_string(),
            row.sample.clone(),
            optional(row.group.clone()),
            row.gene.clone(),
            row.isoform_id.clone(),
            row.site_id.clone(),
            row.site.chrom.clone(),
            row.site.pos0.to_string(),
            row.site.strand.as_char().to_string(),
            row.site.mod_code.clone(),
            optional(row.context.clone()),
            row.site_state.to_string(),
            row.coverage_basis.to_string(),
            row.n_assigned.to_string(),
            optional(row.n_covering),
            row.n_candidate.to_string(),
            row.n_callable.to_string(),
            row.n_modified.to_string(),
            row.n_unmodified.to_string(),
            row.n_unknown.to_string(),
            optional(row.mod_fraction),
            optional(row.mean_probability),
            optional(row.ci_low),
            optional(row.ci_high),
            if row.eligibility_reason.is_eligible() {
                "eligible"
            } else {
                "ineligible"
            }
            .to_owned(),
            row.eligibility_reason.to_string(),
        ])?;
    }
    output.flush().context("flush isoform modification sites")?;
    Ok(())
}

/// Write compact comparison-ready integer counts as deterministic TSV.
pub fn write_isoform_mod_design_tsv<W: Write>(
    writer: W,
    rows: &[IsoformModDesignRow],
) -> anyhow::Result<()> {
    let mut output = WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_writer(writer);
    output.write_record([
        "assay_id",
        "analysis_threshold",
        "sample",
        "group",
        "gene",
        "site_id",
        "mod_code",
        "isoform_id",
        "n_modified",
        "n_unmodified",
        "mod_fraction",
        "eligibility",
        "eligibility_reason",
    ])?;
    for row in rows {
        output.write_record([
            row.assay_id.clone(),
            row.analysis_threshold.to_string(),
            row.sample.clone(),
            optional(row.group.clone()),
            row.gene.clone(),
            row.site_id.clone(),
            row.mod_code.clone(),
            row.isoform_id.clone(),
            row.n_modified.to_string(),
            row.n_unmodified.to_string(),
            optional(row.mod_fraction),
            if row.eligibility_reason.is_eligible() {
                "eligible"
            } else {
                "ineligible"
            }
            .to_owned(),
            row.eligibility_reason.to_string(),
        ])?;
    }
    output
        .flush()
        .context("flush isoform modification design")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::io::coverage::{BamCoverageResult, ReadCoverage};
    use crate::io::mod_calls::ObservationReadResult;
    use crate::model::{Bed12Attrs, Interval};
    use crate::modification::{ModObservation, ModObservationKey, MODIFICATION_SCHEMA_VERSION};

    use super::*;

    fn transcript(name: &str, gene: &str, exons: &[(u32, u32)]) -> Transcript {
        let intervals = exons
            .iter()
            .map(|&(start, end)| Interval::new(Coord::new(start), Coord::new(end)).unwrap())
            .collect::<Vec<_>>();
        Transcript::new(
            "chr1".to_owned(),
            Strand::Plus,
            intervals[0].start,
            intervals.last().unwrap().end,
            name.to_owned(),
            intervals,
            Bed12Attrs {
                score: 0,
                thick_start: Coord::new(exons[0].0),
                thick_end: Coord::new(exons.last().unwrap().1),
                item_rgb: "0".to_owned(),
                extra_fields: vec![
                    "none".to_owned(),
                    "none".to_owned(),
                    "none".to_owned(),
                    "none".to_owned(),
                    "isoform".to_owned(),
                    gene.to_owned(),
                    "none".to_owned(),
                    "none".to_owned(),
                ],
            },
        )
        .unwrap()
    }

    fn metadata(complete: bool) -> AssayMetadata {
        AssayMetadata {
            schema_version: MODIFICATION_SCHEMA_VERSION,
            assay_id: "a1".to_owned(),
            caller: "test".to_owned(),
            caller_version: "1".to_owned(),
            model_id: "model".to_owned(),
            chemistry: "RNA004".to_owned(),
            candidate_rule: "DRACH".to_owned(),
            source_emission_threshold: Some(0.1),
            source_site_filter: "none".to_owned(),
            candidate_observations_complete: complete,
            implicit_skip_policy: ImplicitSkipPolicy::LowProbability,
            coordinate_source: "genome".to_owned(),
            read_id_mapping: "test".to_owned(),
            source_files: Vec::new(),
        }
    }

    fn observation(
        read: &str,
        pos0: u32,
        state: ObservationState,
        probability: Option<f64>,
    ) -> ModObservation {
        ModObservation {
            key: ModObservationKey {
                assay_id: "a1".to_owned(),
                sample: "S1".to_owned(),
                read_id: read.to_owned(),
                site: ModSiteKey {
                    chrom: "chr1".to_owned(),
                    pos0,
                    strand: Strand::Plus,
                    mod_code: "A+a".to_owned(),
                },
            },
            probability,
            observation_state: state,
            context: Some("DRACH".to_owned()),
            source_transcript_id: None,
            source_pos0: None,
        }
    }

    fn run(complete: bool) -> AggregateResult {
        let samples = vec![SampleRow {
            sample: "S1".to_owned(),
            group: Some("control".to_owned()),
            reads: "unused.bed".into(),
        }];
        let isoforms = vec![
            transcript("iso1", "GENE1", &[(100, 120), (200, 220)]),
            transcript("iso2", "GENE1", &[(100, 120), (300, 320)]),
        ];
        let assignments = vec![
            ("S1::r1".to_owned(), "iso1".to_owned()),
            ("S1::r2".to_owned(), "iso1".to_owned()),
            ("S1::r3".to_owned(), "iso1".to_owned()),
        ];
        let observations = vec![
            observation(
                "S1::r1",
                110,
                ObservationState::ExplicitProbability,
                Some(0.9),
            ),
            observation(
                "S1::r2",
                110,
                ObservationState::ExplicitProbability,
                Some(0.2),
            ),
            observation("S1::r3", 110, ObservationState::Unknown, None),
            observation(
                "S1::r1",
                210,
                ObservationState::ImplicitBelowEmissionThreshold,
                None,
            ),
        ];
        let inputs = vec![ModSampleInput {
            sample: "S1".to_owned(),
            assay_id: "a1".to_owned(),
            metadata: metadata(complete),
            observations: ObservationReadResult {
                input_rows: observations.len(),
                duplicate_exact: 0,
                observations,
            },
            coverage: None,
        }];
        aggregate_modifications(
            &samples,
            &isoforms,
            &assignments,
            &inputs,
            &AggregateOptions {
                analysis_thresholds: BTreeMap::from([("a1".to_owned(), 0.5)]),
                min_callable: 1,
                min_read_join_rate: 1.0,
                allow_low_join: false,
                reference_bases: None,
            },
        )
        .unwrap()
    }

    #[test]
    fn exact_bam_coverage_counts_assigned_reads_independently_of_candidates() {
        let samples = vec![SampleRow {
            sample: "S1".to_owned(),
            group: None,
            reads: "unused.bed".into(),
        }];
        let isoforms = vec![transcript("iso1", "GENE1", &[(100, 120)])];
        let assignments = vec![
            ("S1::r1".to_owned(), "iso1".to_owned()),
            ("S1::r2".to_owned(), "iso1".to_owned()),
        ];
        let coverage = BamCoverageResult {
            reads: vec![
                ReadCoverage {
                    read_id: "S1::r1".to_owned(),
                    chrom: "chr1".to_owned(),
                    strand: Strand::Plus,
                    match_blocks: vec![Interval::new(Coord::new(100), Coord::new(120)).unwrap()],
                },
                ReadCoverage {
                    read_id: "S1::r2".to_owned(),
                    chrom: "chr1".to_owned(),
                    strand: Strand::Plus,
                    match_blocks: vec![Interval::new(Coord::new(105), Coord::new(115)).unwrap()],
                },
            ],
            total_records: 2,
            ..BamCoverageResult::default()
        };
        let inputs = vec![ModSampleInput {
            sample: "S1".to_owned(),
            assay_id: "a1".to_owned(),
            metadata: metadata(true),
            observations: ObservationReadResult {
                input_rows: 1,
                duplicate_exact: 0,
                observations: vec![observation(
                    "S1::r1",
                    110,
                    ObservationState::ExplicitProbability,
                    Some(0.9),
                )],
            },
            coverage: Some(coverage),
        }];
        let result = aggregate_modifications(
            &samples,
            &isoforms,
            &assignments,
            &inputs,
            &AggregateOptions {
                analysis_thresholds: BTreeMap::from([("a1".to_owned(), 0.5)]),
                min_callable: 1,
                min_read_join_rate: 1.0,
                allow_low_join: false,
                reference_bases: None,
            },
        )
        .unwrap();
        assert_eq!(result.sites.len(), 1);
        assert_eq!(result.sites[0].coverage_basis, CoverageBasis::BamExact);
        assert_eq!(result.sites[0].n_assigned, 2);
        assert_eq!(result.sites[0].n_covering, Some(2));
        assert_eq!(result.sites[0].n_candidate, 1);
    }

    #[test]
    fn wilson_interval_contains_boundary_fractions_exactly() {
        assert_eq!(wilson_interval(0, 3).unwrap().0, 0.0);
        assert_eq!(wilson_interval(3, 3).unwrap().1, 1.0);
    }

    #[test]
    fn aggregation_distinguishes_zero_na_unknown_and_structural_absence() {
        let result = run(true);
        let shared_iso1 = result
            .sites
            .iter()
            .find(|row| row.isoform_id == "iso1" && row.site.pos0 == 110)
            .unwrap();
        assert_eq!(
            (
                shared_iso1.n_candidate,
                shared_iso1.n_callable,
                shared_iso1.n_modified,
                shared_iso1.n_unmodified,
                shared_iso1.n_unknown
            ),
            (3, 2, 1, 1, 1)
        );
        assert_eq!(shared_iso1.mod_fraction, None);
        assert_eq!(
            shared_iso1.eligibility_reason,
            EligibilityReason::UnknownDenominator
        );

        let shared_iso2 = result
            .sites
            .iter()
            .find(|row| row.isoform_id == "iso2" && row.site.pos0 == 110)
            .unwrap();
        assert_eq!(shared_iso2.n_callable, 0);
        assert_eq!(shared_iso2.mod_fraction, None);
        assert_eq!(
            shared_iso2.eligibility_reason,
            EligibilityReason::LowCallable
        );

        let specific_iso1 = result
            .sites
            .iter()
            .find(|row| row.isoform_id == "iso1" && row.site.pos0 == 210)
            .unwrap();
        assert_eq!(specific_iso1.n_callable, 1);
        assert_eq!(specific_iso1.n_modified, 0);
        assert_eq!(specific_iso1.mod_fraction, Some(0.0));

        let absent_iso2 = result
            .sites
            .iter()
            .find(|row| row.isoform_id == "iso2" && row.site.pos0 == 210)
            .unwrap();
        assert_eq!(absent_iso2.site_state, SiteState::StructurallyAbsent);
        assert_eq!(absent_iso2.mod_fraction, None);
        assert_eq!(
            absent_iso2.eligibility_reason,
            EligibilityReason::SiteAbsent
        );
    }

    #[test]
    fn junction_adjacent_context_is_descriptive_but_ineligible() {
        let samples = vec![SampleRow {
            sample: "S1".to_owned(),
            group: None,
            reads: "unused".into(),
        }];
        let isoforms = vec![transcript("iso1", "GENE1", &[(100, 120), (200, 220)])];
        let assignments = vec![("S1::r1".to_owned(), "iso1".to_owned())];
        let inputs = vec![ModSampleInput {
            sample: "S1".to_owned(),
            assay_id: "a1".to_owned(),
            metadata: metadata(true),
            observations: ObservationReadResult {
                input_rows: 1,
                duplicate_exact: 0,
                observations: vec![observation(
                    "S1::r1",
                    100,
                    ObservationState::ExplicitProbability,
                    Some(0.9),
                )],
            },
            coverage: None,
        }];
        let result = aggregate_modifications(
            &samples,
            &isoforms,
            &assignments,
            &inputs,
            &AggregateOptions {
                analysis_thresholds: BTreeMap::from([("a1".to_owned(), 0.5)]),
                min_callable: 1,
                min_read_join_rate: 1.0,
                allow_low_join: false,
                reference_bases: None,
            },
        )
        .unwrap();

        let row = &result.sites[0];
        assert_eq!(row.site_state, SiteState::ContextDependent);
        assert_eq!((row.n_callable, row.n_modified), (1, 1));
        assert_eq!(row.mod_fraction, Some(1.0));
        assert_eq!(row.eligibility_reason, EligibilityReason::ContextDependent);
    }

    #[test]
    fn low_join_fails_by_default_or_emits_ineligible_rows_when_allowed() {
        let samples = vec![SampleRow {
            sample: "S1".to_owned(),
            group: None,
            reads: "unused".into(),
        }];
        let isoforms = vec![transcript("iso1", "GENE1", &[(100, 120)])];
        let assignments = vec![("S1::r1".to_owned(), "iso1".to_owned())];
        let observations = vec![
            observation(
                "S1::r1",
                110,
                ObservationState::ExplicitProbability,
                Some(0.9),
            ),
            observation(
                "S1::missing",
                110,
                ObservationState::ExplicitProbability,
                Some(0.1),
            ),
        ];
        let inputs = vec![ModSampleInput {
            sample: "S1".to_owned(),
            assay_id: "a1".to_owned(),
            metadata: metadata(true),
            observations: ObservationReadResult {
                input_rows: observations.len(),
                duplicate_exact: 0,
                observations,
            },
            coverage: None,
        }];
        let mut options = AggregateOptions {
            analysis_thresholds: BTreeMap::from([("a1".to_owned(), 0.5)]),
            min_callable: 1,
            min_read_join_rate: 0.75,
            allow_low_join: false,
            reference_bases: None,
        };

        let error = aggregate_modifications(&samples, &isoforms, &assignments, &inputs, &options)
            .unwrap_err()
            .to_string();
        assert!(error.contains("read join rate 0.500000"), "{error}");

        options.allow_low_join = true;
        let result =
            aggregate_modifications(&samples, &isoforms, &assignments, &inputs, &options).unwrap();
        assert_eq!(result.join_qc[0].read_join_rate, 0.5);
        assert_eq!(result.join_qc[0].unknown_read, 1);
        assert_eq!(
            result.sites[0].eligibility_reason,
            EligibilityReason::JoinRateLow
        );
    }

    #[test]
    fn site_local_low_join_is_ineligible_even_when_global_join_passes() {
        let samples = vec![SampleRow {
            sample: "S1".to_owned(),
            group: None,
            reads: "unused".into(),
        }];
        let isoforms = vec![transcript("iso1", "GENE1", &[(100, 320)])];
        let mut assignments = Vec::new();
        let mut observations = Vec::new();
        for index in 0..100 {
            let read = format!("S1::bulk-{index}");
            assignments.push((read.clone(), "iso1".to_owned()));
            observations.push(observation(
                &read,
                110,
                ObservationState::ExplicitProbability,
                Some(0.9),
            ));
        }
        assignments.push(("S1::local-joined".to_owned(), "iso1".to_owned()));
        observations.push(observation(
            "S1::local-joined",
            210,
            ObservationState::ExplicitProbability,
            Some(0.9),
        ));
        for index in 0..3 {
            observations.push(observation(
                &format!("S1::local-missing-{index}"),
                210,
                ObservationState::ExplicitProbability,
                Some(0.1),
            ));
        }
        observations.push(observation(
            "S1::zero-join",
            310,
            ObservationState::ExplicitProbability,
            Some(0.1),
        ));
        let inputs = vec![ModSampleInput {
            sample: "S1".to_owned(),
            assay_id: "a1".to_owned(),
            metadata: metadata(true),
            observations: ObservationReadResult {
                input_rows: observations.len(),
                duplicate_exact: 0,
                observations,
            },
            coverage: None,
        }];
        let result = aggregate_modifications(
            &samples,
            &isoforms,
            &assignments,
            &inputs,
            &AggregateOptions {
                analysis_thresholds: BTreeMap::from([("a1".to_owned(), 0.5)]),
                min_callable: 1,
                min_read_join_rate: 0.9,
                allow_low_join: false,
                reference_bases: None,
            },
        )
        .unwrap();

        assert!(result.join_qc[0].read_join_rate > 0.9);
        let local = result
            .sites
            .iter()
            .find(|row| row.site.pos0 == 210)
            .unwrap();
        assert_eq!(local.eligibility_reason, EligibilityReason::SiteJoinRateLow);
        let zero_join = result
            .site_join_qc
            .iter()
            .find(|row| row.site.pos0 == 310)
            .unwrap();
        assert_eq!(zero_join.joined_rows, 0);
        assert!(!zero_join.passes_min_join_rate);
    }

    #[test]
    fn reference_base_mismatch_is_audited_and_not_callable() {
        let samples = vec![SampleRow {
            sample: "S1".to_owned(),
            group: None,
            reads: "unused".into(),
        }];
        let isoforms = vec![transcript("iso1", "GENE1", &[(100, 120)])];
        let assignments = vec![("S1::r1".to_owned(), "iso1".to_owned())];
        let input_observation = observation(
            "S1::r1",
            110,
            ObservationState::ExplicitProbability,
            Some(0.9),
        );
        let site = input_observation.key.site.clone();
        let inputs = vec![ModSampleInput {
            sample: "S1".to_owned(),
            assay_id: "a1".to_owned(),
            metadata: metadata(true),
            observations: ObservationReadResult {
                input_rows: 1,
                duplicate_exact: 0,
                observations: vec![input_observation],
            },
            coverage: None,
        }];
        let result = aggregate_modifications(
            &samples,
            &isoforms,
            &assignments,
            &inputs,
            &AggregateOptions {
                analysis_thresholds: BTreeMap::from([("a1".to_owned(), 0.5)]),
                min_callable: 1,
                min_read_join_rate: 1.0,
                allow_low_join: false,
                reference_bases: Some(BTreeMap::from([(site, b'C')])),
            },
        )
        .unwrap();

        let row = &result.sites[0];
        assert_eq!(row.site_state, SiteState::ReferenceBaseMismatch);
        assert_eq!((row.n_candidate, row.n_callable, row.n_unknown), (1, 0, 1));
        assert_eq!(row.mod_fraction, None);
        assert_eq!(
            row.eligibility_reason,
            EligibilityReason::ReferenceBaseMismatch
        );
    }

    #[test]
    fn incomplete_candidate_universe_suppresses_fraction_but_keeps_counts_and_scores() {
        let result = run(false);
        let row = result
            .sites
            .iter()
            .find(|row| row.isoform_id == "iso1" && row.site.pos0 == 110)
            .unwrap();
        assert_eq!((row.n_modified, row.n_unmodified), (1, 1));
        assert_eq!(row.mean_probability, Some(0.55));
        assert_eq!(row.mod_fraction, None);
        assert_eq!(
            row.eligibility_reason,
            EligibilityReason::IncompleteCandidateUniverse
        );
    }

    #[test]
    fn sample_level_qc_does_not_split_an_otherwise_compatible_assay() {
        let complete = metadata(true);
        let mut sample_specific = metadata(false);
        sample_specific.implicit_skip_policy = ImplicitSkipPolicy::Unknown;
        assert!(assay_compatible(&complete, &sample_specific));

        let mut different_model = sample_specific;
        different_model.model_id = "other-model".to_owned();
        assert!(!assay_compatible(&complete, &different_model));
    }

    #[test]
    fn rejects_nonunique_mapping() {
        let samples = vec![SampleRow {
            sample: "S1".to_owned(),
            group: None,
            reads: "unused".into(),
        }];
        let isoforms = vec![
            transcript("iso1", "GENE1", &[(100, 120)]),
            transcript("iso2", "GENE1", &[(100, 120)]),
        ];
        let duplicate = vec![
            ("S1::r1".to_owned(), "iso1".to_owned()),
            ("S1::r1".to_owned(), "iso2".to_owned()),
        ];
        let options = AggregateOptions {
            analysis_thresholds: BTreeMap::from([("a1".to_owned(), 0.5)]),
            ..AggregateOptions::default()
        };
        let inputs = vec![ModSampleInput {
            sample: "S1".to_owned(),
            assay_id: "a1".to_owned(),
            metadata: metadata(true),
            observations: ObservationReadResult {
                input_rows: 1,
                duplicate_exact: 0,
                observations: vec![observation(
                    "S1::r1",
                    110,
                    ObservationState::ExplicitProbability,
                    Some(0.9),
                )],
            },
            coverage: None,
        }];
        let error = aggregate_modifications(&samples, &isoforms, &duplicate, &inputs, &options)
            .unwrap_err()
            .to_string();
        assert!(error.contains("maps to multiple isoforms"), "{error}");
    }

    #[test]
    fn rejects_analysis_threshold_below_source_emission_threshold() {
        let samples = vec![SampleRow {
            sample: "S1".to_owned(),
            group: None,
            reads: "unused".into(),
        }];
        let isoforms = vec![transcript("iso1", "GENE1", &[(100, 120)])];
        let assignments = vec![("S1::r1".to_owned(), "iso1".to_owned())];
        let inputs = vec![ModSampleInput {
            sample: "S1".to_owned(),
            assay_id: "a1".to_owned(),
            metadata: metadata(true),
            observations: ObservationReadResult {
                input_rows: 1,
                duplicate_exact: 0,
                observations: vec![observation(
                    "S1::r1",
                    110,
                    ObservationState::ExplicitProbability,
                    Some(0.9),
                )],
            },
            coverage: None,
        }];
        let options = AggregateOptions {
            analysis_thresholds: BTreeMap::from([("a1".to_owned(), 0.05)]),
            ..AggregateOptions::default()
        };
        let error = aggregate_modifications(&samples, &isoforms, &assignments, &inputs, &options)
            .unwrap_err()
            .to_string();
        assert!(
            error.contains("below source emission threshold 0.1"),
            "{error}"
        );
    }

    #[test]
    fn writers_emit_na_and_stable_headers() {
        let result = run(true);
        let mut sites = Vec::new();
        write_isoform_mod_sites_tsv(&mut sites, &result.sites).unwrap();
        let sites = String::from_utf8(sites).unwrap();
        assert!(sites.starts_with("assay_id\tanalysis_threshold\tsample\tgroup\tgene"));
        assert!(sites.contains("\tunavailable\t"));
        assert!(sites.contains("\tNA\t"));

        let mut design = Vec::new();
        write_isoform_mod_design_tsv(&mut design, &result.design).unwrap();
        assert!(String::from_utf8(design)
            .unwrap()
            .contains("\teligible\tok\n"));
    }

    #[test]
    fn generic_centered_iupac_rules_define_structural_flanks() {
        assert_eq!(context_flank(Some("DRACH"), "unused"), 2);
        assert_eq!(context_flank(Some("NNRAYNN"), "unused"), 3);
        assert_eq!(context_flank(Some("NN?NN"), "unused"), 0);
        assert_eq!(context_flank(Some("NNNN"), "unused"), 0);
    }
}
