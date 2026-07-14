use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

mod common;

use common::TestDir;

fn repo_path(relative: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(relative)
}

fn fresh_temp_dir(label: &str) -> TestDir {
    TestDir::new(&format!("output-root-ownership-{label}"))
}

fn assert_rejected_without_clobber(output: &Output, context: &str) {
    assert!(
        !output.status.success(),
        "{context} unexpectedly succeeded\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("pipeline-owned output_root"),
        "{context} did not report the managed-tree policy\nstderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

fn run_initial_flow(output_root: &Path, prefix: &str) {
    let output = Command::new(env!("CARGO_BIN_EXE_trackcluster"))
        .args([
            "flow",
            "--reads",
            repo_path("tests/fixtures/reads.bed").to_str().unwrap(),
            "--reference",
            repo_path("tests/fixtures/ref.bed").to_str().unwrap(),
            "--output-root",
            output_root.to_str().unwrap(),
            "--prefix",
            prefix,
            "--threads",
            "1",
            "--max-reads-per-gene",
            "0",
            "--assignment-mode",
            "fractional",
            "--force",
        ])
        .output()
        .expect("run initial flow");
    assert!(
        output.status.success(),
        "initial flow failed\nstderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn preparedir_rejects_a_reads_input_beneath_output_root() {
    let root = fresh_temp_dir("preparedir");
    let reads = root.join("sample_novel.bed");
    fs::copy(repo_path("tests/fixtures/reads.bed"), &reads).expect("copy reads");
    let original = fs::read(&reads).expect("read original reads");

    let output = Command::new(env!("CARGO_BIN_EXE_trackcluster"))
        .args([
            "preparedir",
            "--reads",
            reads.to_str().unwrap(),
            "--reference",
            repo_path("tests/fixtures/ref.bed").to_str().unwrap(),
            "--output-root",
            root.to_str().unwrap(),
            "--prefix",
            "sample",
        ])
        .output()
        .expect("run preparedir alias regression");

    assert_rejected_without_clobber(&output, "preparedir reads alias");
    assert_eq!(fs::read(&reads).unwrap(), original);
    assert!(!root.join("sample_dedup.bed").exists());
}

#[test]
fn flow_rejects_a_reference_beneath_output_root() {
    let root = fresh_temp_dir("flow-reference");
    let reference = root.join("sample_unused.bed");
    fs::copy(repo_path("tests/fixtures/ref.bed"), &reference).expect("copy reference");
    let original = fs::read(&reference).expect("read original reference");

    let output = Command::new(env!("CARGO_BIN_EXE_trackcluster"))
        .args([
            "flow",
            "--reads",
            repo_path("tests/fixtures/reads.bed").to_str().unwrap(),
            "--reference",
            reference.to_str().unwrap(),
            "--output-root",
            root.to_str().unwrap(),
            "--prefix",
            "sample",
            "--max-reads-per-gene",
            "0",
        ])
        .output()
        .expect("run flow reference alias regression");

    assert_rejected_without_clobber(&output, "flow reference alias");
    assert_eq!(fs::read(&reference).unwrap(), original);
    assert!(!root.join("sample_gene.txt").exists());
}

#[test]
fn flow_rejects_a_manifest_beneath_output_root() {
    let root = fresh_temp_dir("flow-manifest");
    let manifest = root.join("sample_gene.txt");
    fs::copy(repo_path("tests/fixtures/samples.tsv"), &manifest).expect("copy manifest");
    fs::copy(
        repo_path("tests/fixtures/S1.reads.bed"),
        root.join("S1.reads.bed"),
    )
    .expect("copy S1 reads");
    fs::copy(
        repo_path("tests/fixtures/S2.reads.bed"),
        root.join("S2.reads.bed"),
    )
    .expect("copy S2 reads");
    let original = fs::read(&manifest).expect("read original manifest");

    let output = Command::new(env!("CARGO_BIN_EXE_trackcluster"))
        .args([
            "flow",
            "--manifest",
            manifest.to_str().unwrap(),
            "--reference",
            repo_path("tests/fixtures/ref.bed").to_str().unwrap(),
            "--output-root",
            root.to_str().unwrap(),
            "--prefix",
            "sample",
            "--max-reads-per-gene",
            "0",
        ])
        .output()
        .expect("run flow manifest alias regression");

    assert_rejected_without_clobber(&output, "flow manifest alias");
    assert_eq!(fs::read(&manifest).unwrap(), original);
}

#[cfg(unix)]
#[test]
fn flow_rejects_a_pooled_output_hard_linked_to_sample_reads() {
    let parent = fresh_temp_dir("flow-pooled-hardlink");
    let input_root = parent.join("input");
    let output_root = parent.join("output");
    fs::create_dir_all(&input_root).unwrap();
    fs::create_dir_all(&output_root).unwrap();
    let first_reads = input_root.join("S1.reads.bed");
    let second_reads = input_root.join("S2.reads.bed");
    fs::copy(repo_path("tests/fixtures/S1.reads.bed"), &first_reads).unwrap();
    fs::copy(repo_path("tests/fixtures/S2.reads.bed"), &second_reads).unwrap();
    let manifest = input_root.join("manifest.tsv");
    fs::write(
        &manifest,
        format!(
            "sample\tgroup\treads\nS1\tcontrol\t{}\nS2\ttreated\t{}\n",
            first_reads.display(),
            second_reads.display()
        ),
    )
    .unwrap();
    let pooled = output_root.join("sample_pooled_reads.bed");
    fs::hard_link(&first_reads, &pooled).expect("hard-link pooled path to S1 reads");
    let original = fs::read(&first_reads).expect("read original S1 reads");

    let output = Command::new(env!("CARGO_BIN_EXE_trackcluster"))
        .args([
            "flow",
            "--manifest",
            manifest.to_str().unwrap(),
            "--reference",
            repo_path("tests/fixtures/ref.bed").to_str().unwrap(),
            "--output-root",
            output_root.to_str().unwrap(),
            "--prefix",
            "sample",
            "--emit-pooled-reads",
            "--max-reads-per-gene",
            "0",
        ])
        .output()
        .expect("run pooled hard-link regression");

    assert_rejected_without_clobber(&output, "flow pooled/sample hard-link alias");
    assert_eq!(fs::read(&first_reads).unwrap(), original);
    assert_eq!(fs::read(&pooled).unwrap(), original);
    assert!(!output_root.join("sample_gene.txt").exists());
}

#[test]
fn count_output_root_rejects_a_reference_at_a_generated_path() {
    let root = fresh_temp_dir("count-output-root");
    run_initial_flow(&root, "sample");
    let reference = root.join("sample_unused.bed");
    fs::copy(repo_path("tests/fixtures/ref.bed"), &reference).expect("place reference in output");
    let original = fs::read(&reference).expect("read original reference");

    let output = Command::new(env!("CARGO_BIN_EXE_trackcluster"))
        .args([
            "count",
            "--output-root",
            root.to_str().unwrap(),
            "--prefix",
            "sample",
            "--reference",
            reference.to_str().unwrap(),
            "--assignment-mode",
            "fractional",
        ])
        .output()
        .expect("run count output-root alias regression");

    assert_rejected_without_clobber(&output, "count output-root reference alias");
    assert_eq!(fs::read(&reference).unwrap(), original);
}
