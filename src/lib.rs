//! Pure-Rust long-read isoform clustering and counting pipeline.
//!
//! `trackcluster-rs` replaces the Python TrackCluster with native
//! sort/intersect/cluster primitives and no runtime dependency on external CLI tools.
//!
//! # Modules
//!
//! - [`model`] -- core data types (`Transcript`, `Interval`, `Coord`, `Strand`)
//! - [`io`] -- BED12 reading and writing
//! - [`interval`] -- interval sorting, partitioning, intersection, and clustering
//! - [`cluster`] -- isoform clustering algorithms (junction-based and overlap-based)
//! - [`count`] -- expression counting
//! - [`annotate`] -- gene assignment and isoform classification
//! - [`flow`] -- end-to-end pipeline workflows
//!
//! TrackCluster-rs is primarily a command-line application. The supported Rust
//! API is the data-model, BED/manifest I/O, and interval modules. The remaining
//! public modules are compatibility surfaces used by the packaged binaries and
//! may evolve between minor releases.
//!
//! # Reading BED12
//!
//! ```no_run
//! use trackcluster_rs::io::bed::{read_bed12, BedError};
//! use trackcluster_rs::model::Transcript;
//!
//! fn load(path: &std::path::Path) -> Result<Vec<Transcript>, BedError> {
//!     read_bed12(path)?.collect()
//! }
//! ```
//!
//! # Cargo features
//!
//! - `index-binned` enables the reusable fixed-bin interval index. The default
//!   sweep implementation remains available with or without this feature.

#![warn(missing_docs)]

#[doc(hidden)]
#[allow(missing_docs)]
pub mod annotate;
#[allow(missing_docs)]
mod cli;
#[doc(hidden)]
#[allow(missing_docs)]
pub mod cluster;
#[doc(hidden)]
#[allow(missing_docs)]
pub mod config;
#[doc(hidden)]
#[allow(missing_docs)]
pub mod count;
#[doc(hidden)]
#[allow(missing_docs)]
pub mod flow;
mod identity;
/// Safe interval sorting, partitioning, intersection, and locus clustering.
pub mod interval;
/// Strict BED12 and sample-manifest readers and writers.
pub mod io;
mod matching;
/// Validated genomic coordinates, intervals, strands, and transcripts.
pub mod model;
#[allow(missing_docs)]
mod sample;

/// Parse command-line arguments from the process environment and run the CLI.
///
/// Packaged binaries use this narrow entry point so the internal Clap command
/// representation does not become part of the supported Rust API.
pub fn run_cli_from_env() -> anyhow::Result<()> {
    cli::run_from_env()
}
