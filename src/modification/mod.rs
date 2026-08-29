//! Isoform-level RNA modification observations, aggregation, and contrasts.

/// Unique-assignment sample/isoform/site aggregation and TSV outputs.
pub mod aggregate;
/// Explicit shared-site effect-only contrasts.
pub mod contrast;
/// Transactional output generations for flow-integrated modification analysis.
pub(crate) mod generation;
/// Streaming per-site summaries of isoform modification audit tables.
pub mod site_summary;
/// Deterministic molecule-level pseudo-sample bundle generation.
pub mod subsample;
/// Modification data contracts shared by caller importers and aggregation.
pub mod types;

pub use types::{
    AssayMetadata, CoverageBasis, EligibilityProfile, EligibilityReason, ImplicitSkipPolicy,
    ModObservation, ModObservationKey, ModSiteKey, ObservationState, ProvenanceStatus, SiteState,
    MODIFICATION_SCHEMA_VERSION,
};
