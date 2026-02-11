//! Pure-Rust long-read isoform clustering and counting pipeline.
//!
//! `trackcluster-rs` replaces the Python TrackCluster with native
//! sort/intersect/cluster primitives and no runtime `bedtools` dependency.
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
//! - [`cli`] -- command-line interface

pub mod annotate;
pub mod cli;
pub mod cluster;
pub mod count;
pub mod flow;
pub mod interval;
pub mod io;
pub mod model;
