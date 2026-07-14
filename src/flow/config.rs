//! Composable, validated configuration for the end-to-end flow.

use crate::cluster::clusterj::{
    JunctionCorrectionOptions, Name2Mode, PlatformPreset, ResolvedPlatformOptions, SlMergeOptions,
    ThreePrimeMergeOptions,
};
use crate::count::{AssignmentMode, UniqueAssignmentOptions};

/// Overlap fractions used while assigning input records to per-gene folders.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PrepareConfig {
    /// Minimum fraction of a read span that must overlap a reference span.
    pub fraction_read: f64,
    /// Minimum fraction of a reference span that must overlap a read span.
    pub fraction_ref: f64,
}

impl Default for PrepareConfig {
    fn default() -> Self {
        Self {
            fraction_read: 0.01,
            fraction_ref: 0.05,
        }
    }
}

impl PrepareConfig {
    /// Validate both overlap fractions.
    pub fn validate(self) -> Result<(), crate::config::ParameterError> {
        crate::config::UnitFraction::new("prepare read overlap fraction", self.fraction_read)?;
        crate::config::UnitFraction::new("prepare reference overlap fraction", self.fraction_ref)?;
        Ok(())
    }
}

/// Optional CLI overrides applied on top of a platform preset.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct JunctionOverrides {
    /// Internal junction-correction offset.
    pub correction_offset: Option<u32>,
    /// Minimum weighted support for an uncorrected junction.
    pub correction_min_support: Option<u32>,
    /// SL-supported partial-read 5-prime offset.
    pub sl_partial_five_prime_offset: Option<u32>,
    /// SL-supported same-junction 5-prime offset.
    pub sl_same_junction_five_prime_offset: Option<u32>,
    /// SL 5-prime cluster offset.
    pub sl_five_prime_cluster_offset: Option<u32>,
    /// Minimum support for an SL 5-prime cluster.
    pub sl_five_prime_min_support: Option<usize>,
    /// Same-junction 3-prime offset.
    pub same_junction_three_prime_offset: Option<u32>,
    /// 3-prime cluster offset.
    pub three_prime_cluster_offset: Option<u32>,
    /// Minimum support for a 3-prime cluster.
    pub three_prime_min_support: Option<usize>,
}

/// Resolved junction-mode scientific settings.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct JunctionConfig {
    /// Preset from which the resolved settings originated.
    pub platform_preset: PlatformPreset,
    /// Internal junction-correction settings.
    pub correction: JunctionCorrectionOptions,
    /// SL/5-prime merge settings.
    pub sl: SlMergeOptions,
    /// 3-prime merge settings.
    pub three_prime: ThreePrimeMergeOptions,
}

impl Default for JunctionConfig {
    fn default() -> Self {
        Self::resolve(PlatformPreset::Generic, JunctionOverrides::default())
    }
}

impl JunctionConfig {
    /// Resolve one platform preset plus explicit CLI overrides into a complete configuration.
    pub fn resolve(preset: PlatformPreset, overrides: JunctionOverrides) -> Self {
        let resolved = crate::cluster::clusterj::resolve_platform_options(
            preset,
            overrides.correction_offset,
            overrides.correction_min_support,
            overrides.sl_partial_five_prime_offset,
            overrides.sl_same_junction_five_prime_offset,
            overrides.sl_five_prime_cluster_offset,
            overrides.sl_five_prime_min_support,
            overrides.same_junction_three_prime_offset,
            overrides.three_prime_cluster_offset,
            overrides.three_prime_min_support,
        );
        Self::from_resolved(preset, resolved)
    }

    /// Pair a preset identity with already-resolved junction settings.
    pub const fn from_resolved(preset: PlatformPreset, resolved: ResolvedPlatformOptions) -> Self {
        Self {
            platform_preset: preset,
            correction: resolved.junction_correction,
            sl: resolved.sl_options,
            three_prime: resolved.three_prime_options,
        }
    }

    /// Return the algorithm-facing resolved settings.
    pub const fn resolved(self) -> ResolvedPlatformOptions {
        ResolvedPlatformOptions {
            junction_correction: self.correction,
            sl_options: self.sl,
            three_prime_options: self.three_prime,
        }
    }

    /// Validate every resolved junction setting.
    pub fn validate(self) -> Result<(), crate::config::ParameterError> {
        self.resolved().validate()
    }
}

/// Scientific settings specific to overlap-mode clustering.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OverlapConfig {
    /// Pass-one distance cutoff.
    pub cutoff1: f64,
    /// Pass-two distance cutoff.
    pub cutoff2: f64,
    /// Relative intron-distance weight.
    pub intron_weight: f64,
}

impl Default for OverlapConfig {
    fn default() -> Self {
        Self {
            cutoff1: crate::cluster::cluster_overlap::DEFAULT_CUTOFF1,
            cutoff2: crate::cluster::cluster_overlap::DEFAULT_CUTOFF2,
            intron_weight: crate::cluster::cluster_overlap::DEFAULT_INTRON_WEIGHT,
        }
    }
}

impl OverlapConfig {
    /// Validate overlap cutoffs and weight.
    pub fn validate(self) -> Result<(), crate::config::ParameterError> {
        crate::config::UnitFraction::new("overlap pass-1 cutoff", self.cutoff1)?;
        crate::config::UnitFraction::new("overlap pass-2 cutoff", self.cutoff2)?;
        crate::config::NonNegativeWeight::new("overlap intron weight", self.intron_weight)?;
        Ok(())
    }
}

/// Settings shared by both clustering algorithms.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ClusteringConfig {
    /// Smith-Waterman score cutoff for SL support.
    pub sw_score: i64,
    /// Maximum records in one bounded merge batch; zero disables batching where supported.
    pub batch_size: usize,
    /// Maximum number of bounded merge rounds.
    pub batch_rounds: usize,
    /// Encoding used for BED `name2` output.
    pub name2_mode: Name2Mode,
    /// Junction-mode configuration.
    pub junction: JunctionConfig,
    /// Overlap-mode configuration.
    pub overlap: OverlapConfig,
}

impl Default for ClusteringConfig {
    fn default() -> Self {
        Self {
            sw_score: crate::cluster::clusterj::DEFAULT_SW_SCORE,
            batch_size: 500,
            batch_rounds: 100,
            name2_mode: Name2Mode::Coverage,
            junction: JunctionConfig::default(),
            overlap: OverlapConfig::default(),
        }
    }
}

impl ClusteringConfig {
    /// Validate algorithm-specific settings before execution begins.
    pub fn validate(self) -> Result<(), crate::config::ParameterError> {
        crate::config::BatchRounds::new(self.batch_rounds)?;
        self.junction.validate()?;
        self.overlap.validate()?;
        Ok(())
    }

    /// Convert the overlap settings into the low-level algorithm configuration.
    pub const fn overlap_options(self) -> crate::cluster::cluster_overlap::ClusterOptions {
        crate::cluster::cluster_overlap::ClusterOptions {
            cutoff1: self.overlap.cutoff1,
            cutoff2: self.overlap.cutoff2,
            intron_weight: self.overlap.intron_weight,
            sw_score: self.sw_score,
            name2_mode: self.name2_mode,
            batch_size: self.batch_size,
            batch_rounds: self.batch_rounds,
        }
    }

    /// Run one validated per-gene clustering operation without positional scientific options.
    pub fn cluster_gene(
        self,
        mode: crate::flow::full::ClusterMode,
        reads: &[crate::model::Transcript],
        references: &[crate::model::Transcript],
        threads: usize,
    ) -> Result<crate::cluster::result::ClusterResult, crate::config::ParameterError> {
        self.validate()?;
        crate::config::WorkerThreads::new(threads)?;
        match mode {
            crate::flow::full::ClusterMode::Clusterj => {
                crate::cluster::clusterj::try_clusterj_with_options(
                    reads,
                    Some(references),
                    threads,
                    self.sw_score,
                    self.batch_size,
                    self.batch_rounds,
                    self.name2_mode,
                    self.junction.sl,
                    self.junction.three_prime,
                    self.junction.correction,
                )
            }
            crate::flow::full::ClusterMode::Cluster => {
                crate::cluster::cluster_overlap::try_cluster_with_options(
                    reads,
                    Some(references),
                    threads,
                    self.overlap_options(),
                )
            }
        }
    }
}

/// Counting and unique-assignment settings.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct CountingConfig {
    /// How reads with multiple candidate isoforms contribute to counts.
    pub assignment_mode: AssignmentMode,
    /// Settings used by unique assignment.
    pub unique_assignment: UniqueAssignmentOptions,
}

/// Per-gene read downsampling settings.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DownsampleConfig {
    /// Restrict downsampling to these genes; empty means every gene.
    pub genes: Vec<String>,
    /// Maximum reads retained per selected gene; zero disables downsampling.
    pub max_reads_per_gene: usize,
    /// Base seed mixed with the biological gene identifier.
    pub seed: u64,
}

impl Default for DownsampleConfig {
    fn default() -> Self {
        Self {
            genes: Vec::new(),
            max_reads_per_gene: 50_000,
            seed: 1,
        }
    }
}

impl DownsampleConfig {
    /// Validate all explicitly selected biological gene identifiers.
    pub fn validate(&self) -> anyhow::Result<()> {
        crate::flow::path_key::validate_gene_ids(self.genes.iter().map(String::as_str))
            .map(|_| ())
            .context("invalid --downsample-gene value")
    }

    /// Return whether one gene should be downsampled.
    pub fn selects(&self, gene: &str) -> bool {
        self.max_reads_per_gene > 0
            && (self.genes.is_empty() || self.genes.iter().any(|selected| selected == gene))
    }
}

/// Policy applied to recoverable failures isolated to one gene.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum GeneErrorPolicy {
    /// Record gene-local failures, exclude those genes, and continue with verified results.
    #[default]
    Continue,
    /// Record every failure, then stop before batch-level downstream publication.
    Strict,
}

/// Policy applied to malformed records in read BED inputs.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum InvalidReadPolicy {
    /// Exclude only the malformed read track and continue processing the gene.
    #[default]
    Skip,
    /// Preserve strict historical behavior and fail the enclosing stage.
    Fail,
}

impl std::fmt::Display for InvalidReadPolicy {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter.write_str(match self {
            Self::Skip => "skip",
            Self::Fail => "fail",
        })
    }
}

impl std::str::FromStr for InvalidReadPolicy {
    type Err = String;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        if value.eq_ignore_ascii_case("skip") {
            return Ok(Self::Skip);
        }
        if value.eq_ignore_ascii_case("fail") {
            return Ok(Self::Fail);
        }
        Err(format!(
            "invalid read policy {value:?}; expected one of: skip, fail"
        ))
    }
}

impl std::fmt::Display for GeneErrorPolicy {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter.write_str(match self {
            Self::Continue => "continue",
            Self::Strict => "strict",
        })
    }
}

/// Bounded execution and progress-reporting settings.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RuntimeConfig {
    /// Maximum number of per-gene worker threads.
    pub threads: usize,
    /// Rebuild artifacts even when their completion manifest is reusable.
    pub force: bool,
    /// Emit one progress line after this many completed genes.
    pub progress_every: usize,
    /// Heartbeat interval in seconds; zero disables heartbeats.
    pub heartbeat_seconds: u64,
    /// Maximum number of in-flight genes shown by a stalled heartbeat.
    pub heartbeat_top: usize,
    /// Policy applied after all gene-local tasks have completed.
    pub gene_error_policy: GeneErrorPolicy,
    /// Policy applied to record-local read parsing and empty read identifiers.
    pub invalid_read_policy: InvalidReadPolicy,
}

impl Default for RuntimeConfig {
    fn default() -> Self {
        Self {
            threads: 1,
            force: false,
            progress_every: 1_000,
            heartbeat_seconds: 60,
            heartbeat_top: 5,
            gene_error_policy: GeneErrorPolicy::Continue,
            invalid_read_policy: InvalidReadPolicy::Skip,
        }
    }
}

impl RuntimeConfig {
    /// Validate runtime bounds before any worker is started.
    pub fn validate(self) -> Result<(), crate::config::ParameterError> {
        crate::config::WorkerThreads::new(self.threads)?;
        Ok(())
    }

    /// Return the number of workers to start for a finite work set.
    pub fn worker_count(self, work_items: usize) -> usize {
        self.threads.min(work_items)
    }
}

use anyhow::Context as _;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn composed_defaults_validate() {
        PrepareConfig::default().validate().unwrap();
        ClusteringConfig::default().validate().unwrap();
        RuntimeConfig::default().validate().unwrap();
        DownsampleConfig::default().validate().unwrap();
    }

    #[test]
    fn junction_overrides_are_resolved_once() {
        let config = JunctionConfig::resolve(
            PlatformPreset::Rna004,
            JunctionOverrides {
                correction_offset: Some(17),
                three_prime_cluster_offset: None,
                ..JunctionOverrides::default()
            },
        );
        assert_eq!(config.correction.offset, 17);
        assert_eq!(config.three_prime.three_prime_cluster_offset, 17);
        config.validate().unwrap();
    }

    #[test]
    fn runtime_worker_count_is_bounded_by_work() {
        let runtime = RuntimeConfig {
            threads: 8,
            ..RuntimeConfig::default()
        };
        assert_eq!(runtime.worker_count(3), 3);
    }

    #[test]
    fn clustering_rejects_zero_batch_rounds() {
        let config = ClusteringConfig {
            batch_rounds: 0,
            ..ClusteringConfig::default()
        };
        let error = config.validate().expect_err("zero rounds must be rejected");
        assert!(error.to_string().contains("batch rounds"));
    }
}
