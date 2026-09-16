mod common;

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

use common::{assert_success, TestDir};

fn repo_path(path: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(path)
}

fn fresh_temp_dir(prefix: &str) -> TestDir {
    TestDir::new(prefix)
}

fn normalized_lines(path: &Path) -> Vec<String> {
    let content = fs::read_to_string(path).unwrap_or_default();
    let lines: Vec<String> = content
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty())
        .map(ToOwned::to_owned)
        .collect();
    lines
}

#[test]
fn overlap_cli_retains_sl_structure_through_ordinary_intermediates() {
    let root = fresh_temp_dir("overlap_sl_intermediate");
    let reads = root.join("reads.bed");
    let reference = root.join("reference.bed");
    fs::write(
        &reference,
        "chr1\t0\t1200\tref\t100\t+\t0\t0\t0\t2\t200,700,\t0,500,\n",
    )
    .unwrap();
    fs::write(
        &reads,
        concat!(
            "chr1\t100\t800\tsl\t20\t+\t0\t0\t0\t2\t100,300,\t0,400,\n",
            "chr1\t100\t1000\tlong_non_sl\t0\t+\t0\t0\t0\t2\t100,500,\t0,400,\n",
            "chr1\t150\t1100\tlate_non_sl\t0\t+\t0\t0\t0\t2\t50,600,\t0,350,\n",
        ),
    )
    .unwrap();
    let mut expected_bed = None;
    for batch in [None, Some(0), Some(1), Some(2), Some(3)] {
        let out = root.join(format!("batch_{batch:?}.bed"));
        let mut command = Command::new(env!("CARGO_BIN_EXE_trackcluster"));
        command
            .args(["cluster", "--reads"])
            .arg(&reads)
            .arg("--reference")
            .arg(&reference)
            .arg("--out")
            .arg(&out);
        if let Some(batch) = batch {
            command.args([
                "--batch-size",
                &batch.to_string(),
                "--threads",
                "2",
                "--sw-score",
                "11",
                "--sl-partial-5prime-offset",
                "15",
            ]);
        }
        assert_success(&command.output().unwrap(), "overlap SL intermediate");
        let isoforms = trackcluster_rs::io::bed::read_bed12(&out)
            .unwrap()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert_eq!(isoforms.len(), 2, "batch={batch:?}");
        let sl = isoforms
            .iter()
            .find(|tx| tx.tx_start.get() == 100 && tx.tx_end.get() == 800)
            .unwrap();
        assert_eq!((sl.tx_start.get(), sl.tx_end.get()), (100, 800));
        assert_eq!(sl.score, 20);
        let mut mapping = normalized_lines(&out.with_extension("read_to_isoform.tsv"));
        mapping.sort();
        assert_eq!(
            mapping,
            [
                "late_non_sl\tref".to_owned(),
                "long_non_sl\tref".to_owned(),
                format!("sl\t{}", sl.name)
            ],
            "batch={batch:?}"
        );
        assert!(normalized_lines(&out.with_extension("unused.bed")).is_empty());
        let mut bed = normalized_lines(&out);
        bed.sort();
        assert_eq!(expected_bed.get_or_insert_with(|| bed.clone()), &bed);
    }
}

#[test]
fn overlap_cli_applies_configured_sl_five_prime_window() {
    let root = fresh_temp_dir("overlap_sl_window");
    let reads = root.join("reads.bed");
    let reference = root.join("reference.bed");
    fs::write(
        &reference,
        "chr1\t90\t101\tref\t100\t+\t0\t0\t0\t1\t11,\t0,\n",
    )
    .unwrap();
    fs::write(
        &reads,
        concat!(
            "chr1\t100\t1000\tlong\t20\t+\t0\t0\t0\t1\t900,\t0,\n",
            "chr1\t115\t900\tshort\t20\t+\t0\t0\t0\t1\t785,\t0,\n",
        ),
    )
    .unwrap();
    for (window, expected_isoforms) in [(14, 3), (15, 2)] {
        let out = root.join(format!("offset_{window}.bed"));
        let output = Command::new(env!("CARGO_BIN_EXE_trackcluster"))
            .args([
                "cluster",
                "--cutoff1",
                "0",
                "--sw-score",
                "11",
                "--sl-partial-5prime-offset",
                &window.to_string(),
            ])
            .arg("-s")
            .arg(&reads)
            .arg("-r")
            .arg(&reference)
            .arg("-o")
            .arg(&out)
            .output()
            .unwrap();
        assert_success(&output, "overlap SL window");
        assert_eq!(normalized_lines(&out).len(), expected_isoforms);
        assert_eq!(
            normalized_lines(&out.with_extension("read_to_isoform.tsv")).len(),
            2
        );
    }
}

#[test]
fn cluster_overlap_plain_bed12_reference_and_unmatched_reads_match_goldens() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");

    let reads = repo_path("tests/fixtures/cluster_overlap/plain_reads.bed");
    let reference = repo_path("tests/fixtures/cluster_overlap/plain_ref.bed");
    let golden_isoforms = repo_path("tests/golden/cluster/plain_isoform.bed");
    let golden_mapping = repo_path("tests/golden/cluster/plain_isoform.read_to_isoform.tsv");
    let golden_unused = repo_path("tests/golden/cluster/plain_isoform.unused.bed");

    let out_dir = fresh_temp_dir("cluster_overlap_plain");
    let out_bed = out_dir.join("isoform.bed");

    let output = Command::new(exe)
        .args([
            "cluster",
            "-s",
            reads.to_str().unwrap(),
            "-r",
            reference.to_str().unwrap(),
            "-o",
            out_bed.to_str().unwrap(),
        ])
        .output()
        .expect("run overlap cluster");
    assert_success(&output, "overlap-cluster golden run");

    let produced_mapping = out_bed.with_extension("read_to_isoform.tsv");
    let produced_unused = out_bed.with_extension("unused.bed");

    assert_eq!(
        normalized_lines(&out_bed),
        normalized_lines(&golden_isoforms)
    );
    assert_eq!(
        normalized_lines(&produced_mapping),
        normalized_lines(&golden_mapping)
    );
    assert_eq!(
        normalized_lines(&produced_unused),
        normalized_lines(&golden_unused)
    );
}
