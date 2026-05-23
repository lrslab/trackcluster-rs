use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

fn repo_path(path: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(path)
}

fn fresh_temp_dir(prefix: &str) -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("time went backwards")
        .as_nanos();
    let mut dir = std::env::temp_dir();
    dir.push(format!(
        "trackcluster_rs_{}_{}_{}",
        prefix,
        std::process::id(),
        nanos
    ));
    fs::create_dir_all(&dir).expect("create temp dir");
    dir
}

fn normalized_lines(path: &Path) -> Vec<String> {
    let content = fs::read_to_string(path).unwrap_or_default();
    let mut lines: Vec<String> = content
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty())
        .map(ToOwned::to_owned)
        .collect();
    lines.sort();
    lines
}

fn count_sum(path: &Path) -> f64 {
    fs::read_to_string(path)
        .expect("read count csv")
        .lines()
        .skip(1)
        .filter(|line| !line.trim().is_empty())
        .map(|line| {
            line.split(',')
                .nth(1)
                .expect("count column")
                .parse::<f64>()
                .expect("parse count")
        })
        .sum()
}

#[test]
fn flow_runs_end_to_end_and_matches_goldens() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");

    let reads = repo_path("tests/fixtures/reads.bed");
    let reference = repo_path("tests/fixtures/ref.bed");
    let golden_isoforms = repo_path("tests/golden/clusterj/isoform.bed");
    let golden_unused = repo_path("tests/golden/clusterj/isoform.unused.bed");
    let golden_count = repo_path("tests/golden/count/isoform_count.csv");

    let out_dir = fresh_temp_dir("flow");
    let prefix = "sample";

    let status = Command::new(exe)
        .args([
            "flow",
            "-s",
            reads.to_str().unwrap(),
            "-r",
            reference.to_str().unwrap(),
            "-o",
            out_dir.to_str().unwrap(),
            "--prefix",
            prefix,
            "--threads",
            "1",
            "--force",
        ])
        .status()
        .expect("run flow");
    assert!(status.success());

    let isoform_out = out_dir.join(format!("{prefix}_isoform.bed"));
    let unused_out = out_dir.join(format!("{prefix}_unused.bed"));
    let count_out = out_dir.join(format!("{prefix}_isoform_count.csv"));
    let unique_mapping_out = out_dir.join(format!("{prefix}_read_to_isoform.unique.tsv"));

    assert_eq!(
        normalized_lines(&isoform_out),
        normalized_lines(&golden_isoforms)
    );
    assert_eq!(
        normalized_lines(&unused_out),
        normalized_lines(&golden_unused)
    );
    assert_eq!(
        normalized_lines(&count_out),
        normalized_lines(&golden_count)
    );
    assert!(unique_mapping_out.exists());
    assert_eq!(count_sum(&count_out), 1.0);

    assert!(out_dir.join("GENEA/GENEA_simple_coveragej.bed").exists());
    assert!(out_dir.join("GENEA/GENEA_unused.bed").exists());
    assert!(out_dir.join("GENEA/GENEA_read_to_isoform.tsv").exists());

    assert!(out_dir.join(format!("{prefix}_desc.txt")).exists());
    assert!(out_dir.join(format!("{prefix}_class4.txt")).exists());
    assert!(out_dir.join(format!("{prefix}_fusion.txt")).exists());
    assert!(out_dir.join(format!("{prefix}_class12.txt")).exists());
}

#[test]
fn flow_count_only_reuses_completed_gene_outputs() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");

    let reads = repo_path("tests/fixtures/reads.bed");
    let reference = repo_path("tests/fixtures/ref.bed");
    let golden_count = repo_path("tests/golden/count/isoform_count.csv");

    let out_dir = fresh_temp_dir("flow_count_only");
    let prefix = "sample";

    let status = Command::new(exe)
        .args([
            "flow",
            "-s",
            reads.to_str().unwrap(),
            "-r",
            reference.to_str().unwrap(),
            "-o",
            out_dir.to_str().unwrap(),
            "--prefix",
            prefix,
            "--threads",
            "1",
            "--force",
        ])
        .status()
        .expect("run initial flow");
    assert!(status.success());

    let count_out = out_dir.join(format!("{prefix}_isoform_count.csv"));
    let unique_mapping_out = out_dir.join(format!("{prefix}_read_to_isoform.unique.tsv"));
    fs::remove_file(&count_out).expect("remove count output");
    fs::remove_file(&unique_mapping_out).expect("remove unique mapping output");

    let status = Command::new(exe)
        .args([
            "flow",
            "--count-only",
            "-r",
            reference.to_str().unwrap(),
            "-o",
            out_dir.to_str().unwrap(),
            "--prefix",
            prefix,
        ])
        .status()
        .expect("run count-only flow");
    assert!(status.success());

    assert_eq!(
        normalized_lines(&count_out),
        normalized_lines(&golden_count)
    );
    assert!(unique_mapping_out.exists());
    assert!(out_dir.join(format!("{prefix}_desc.txt")).exists());
}

#[test]
fn flow_manifest_writes_multi_sample_usage_outputs() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");

    let manifest = repo_path("tests/fixtures/samples.tsv");
    let reference = repo_path("tests/fixtures/ref.bed");

    let out_dir = fresh_temp_dir("flow_manifest");
    let prefix = "pooled";

    let status = Command::new(exe)
        .args([
            "flow",
            "--manifest",
            manifest.to_str().unwrap(),
            "-r",
            reference.to_str().unwrap(),
            "-o",
            out_dir.to_str().unwrap(),
            "--prefix",
            prefix,
            "--emit-pooled-reads",
            "--threads",
            "1",
            "--force",
        ])
        .status()
        .expect("run flow manifest mode");
    assert!(status.success());

    let pooled_reads = out_dir.join(format!("{prefix}_pooled_reads.bed"));
    assert!(pooled_reads.exists());
    let pooled = fs::read_to_string(pooled_reads).expect("read pooled reads");
    assert!(pooled.contains("\tS1::read_s1\t"));
    assert!(pooled.contains("\tS2::read_s2\t"));

    let long_tsv = out_dir.join(format!("{prefix}.isoform_usage.long.tsv"));
    let matrix_tsv = out_dir.join(format!("{prefix}.isoform_counts.matrix.tsv"));
    let multi_count_csv = out_dir.join(format!("{prefix}.isoform_count.csv"));
    let main_count_csv = out_dir.join(format!("{prefix}_isoform_count.csv"));
    let group_tsv = out_dir.join(format!("{prefix}.isoform_usage.group.tsv"));
    assert!(long_tsv.exists());
    assert!(matrix_tsv.exists());
    assert!(multi_count_csv.exists());
    assert!(main_count_csv.exists());
    assert!(group_tsv.exists());

    let long_content = fs::read_to_string(long_tsv).expect("read long tsv");
    let mut sample_prop_sum = std::collections::HashMap::<String, f64>::new();
    for (idx, line) in long_content.lines().enumerate() {
        if idx == 0 {
            assert_eq!(
                line,
                "gene\tisoform_id\tsample\tgroup\tcount\tproportion\tgene_total"
            );
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields.len(), 7);
        let key = format!("{}\t{}", fields[0], fields[2]); // gene + sample
        let proportion: f64 = fields[5].parse().expect("parse proportion");
        *sample_prop_sum.entry(key).or_insert(0.0) += proportion;
    }
    for total in sample_prop_sum.values() {
        assert!((*total - 1.0).abs() < 1e-9);
    }

    let matrix_content = fs::read_to_string(matrix_tsv).expect("read matrix tsv");
    let mut matrix_totals = std::collections::HashMap::<String, f64>::new();
    for (idx, line) in matrix_content.lines().enumerate() {
        if idx == 0 {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let total = fields[2..]
            .iter()
            .map(|value| value.parse::<f64>().expect("parse sample count"))
            .sum();
        matrix_totals.insert(fields[1].to_owned(), total);
    }

    for path in [multi_count_csv, main_count_csv] {
        let csv = fs::read_to_string(path).expect("read aggregate count csv");
        for (idx, line) in csv.lines().enumerate() {
            if idx == 0 {
                assert_eq!(line, "isoform_id,count");
                continue;
            }
            let fields: Vec<&str> = line.split(',').collect();
            assert_eq!(fields.len(), 2);
            let matrix_total = matrix_totals
                .get(fields[0])
                .copied()
                .expect("matrix row for isoform");
            let aggregate = fields[1].parse::<f64>().expect("parse aggregate count");
            assert!((aggregate - matrix_total).abs() < 1e-9);
        }
    }
}

#[test]
fn flow_manifest_skips_pooled_reads_by_default() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");

    let manifest = repo_path("tests/fixtures/samples.tsv");
    let reference = repo_path("tests/fixtures/ref.bed");

    let out_dir = fresh_temp_dir("flow_manifest_no_pool");
    let prefix = "pooled";

    let status = Command::new(exe)
        .args([
            "flow",
            "--manifest",
            manifest.to_str().unwrap(),
            "-r",
            reference.to_str().unwrap(),
            "-o",
            out_dir.to_str().unwrap(),
            "--prefix",
            prefix,
            "--threads",
            "1",
            "--force",
        ])
        .status()
        .expect("run flow manifest mode");
    assert!(status.success());

    let pooled_reads = out_dir.join(format!("{prefix}_pooled_reads.bed"));
    assert!(!pooled_reads.exists());
}

#[test]
fn flow_manifest_downsamples_gene_over_cutoff_and_writes_scale_factors() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");

    let manifest = repo_path("tests/fixtures/samples.tsv");
    let reference = repo_path("tests/fixtures/ref.bed");

    let out_dir = fresh_temp_dir("flow_manifest_downsample");
    let prefix = "pooled";

    let status = Command::new(exe)
        .args([
            "flow",
            "--manifest",
            manifest.to_str().unwrap(),
            "-r",
            reference.to_str().unwrap(),
            "-o",
            out_dir.to_str().unwrap(),
            "--prefix",
            prefix,
            "--threads",
            "1",
            "--force",
            "--max-reads-per-gene",
            "1",
            "--downsample-seed",
            "1",
        ])
        .status()
        .expect("run flow manifest downsample");
    assert!(status.success());

    let summary = out_dir.join("clusterj_batch_downsample.tsv");
    assert!(summary.exists());
    let summary_text = fs::read_to_string(summary).unwrap();
    assert!(summary_text
        .lines()
        .any(|line| line.starts_with("GENEA\t2\t1\t2")));

    assert!(out_dir.join("GENEA/downsample.tsv").exists());

    let main_count = out_dir.join(format!("{prefix}_isoform_count.csv"));
    let multi_count = out_dir.join(format!("{prefix}.isoform_count.csv"));
    assert!((count_sum(&main_count) - 2.0).abs() < 1e-9);
    assert!((count_sum(&multi_count) - 2.0).abs() < 1e-9);
}

#[test]
fn flow_overlap_mode_runs_end_to_end() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");

    let reads = {
        let read1 = fs::read_to_string(repo_path("tests/fixtures/S1.reads.bed")).expect("read S1");
        let read2 = fs::read_to_string(repo_path("tests/fixtures/S2.reads.bed")).expect("read S2");
        let path = fresh_temp_dir("flow_overlap_reads").join("reads.bed");
        fs::write(&path, format!("{read1}{read2}")).expect("write overlap reads");
        path
    };
    let reference = repo_path("tests/fixtures/ref.bed");

    let out_dir = fresh_temp_dir("flow_overlap");
    let prefix = "sample";

    let status = Command::new(exe)
        .args([
            "flow",
            "--cluster-mode",
            "cluster",
            "-s",
            reads.to_str().unwrap(),
            "-r",
            reference.to_str().unwrap(),
            "-o",
            out_dir.to_str().unwrap(),
            "--prefix",
            prefix,
            "--threads",
            "1",
            "--batch-size",
            "1",
            "--batch-rounds",
            "4",
            "--force",
        ])
        .status()
        .expect("run flow overlap mode");
    assert!(status.success());

    let isoform_out = out_dir.join(format!("{prefix}_isoform.bed"));
    let unused_out = out_dir.join(format!("{prefix}_unused.bed"));
    let mapping_out = out_dir.join(format!("{prefix}_read_to_isoform.tsv"));
    let count_out = out_dir.join(format!("{prefix}_isoform_count.csv"));

    assert!(isoform_out.exists());
    assert!(unused_out.exists());
    assert!(mapping_out.exists());
    assert!(count_out.exists());

    let per_gene_unused = out_dir.join("GENEA/GENEA_unused.bed");
    assert!(out_dir.join("GENEA/GENEA_simple_coverage.bed").exists());
    assert!(per_gene_unused.exists());
    assert!(out_dir.join("cluster_batch_summary.txt").exists());

    let mapping = fs::read_to_string(mapping_out).expect("read merged mapping");
    assert!(mapping.lines().any(|line| line.starts_with("read_s1\t")));
    assert!(mapping.lines().any(|line| line.starts_with("read_s2\t")));
    assert!((count_sum(&count_out) - 2.0).abs() < 1e-9);

    assert!(out_dir.join(format!("{prefix}_desc.txt")).exists());
    assert!(out_dir.join(format!("{prefix}_class4.txt")).exists());
    assert!(out_dir.join(format!("{prefix}_fusion.txt")).exists());
    assert!(out_dir.join(format!("{prefix}_class12.txt")).exists());

    fs::remove_file(&per_gene_unused).expect("remove per-gene unused");
    let status = Command::new(exe)
        .args([
            "flow",
            "--cluster-mode",
            "cluster",
            "-s",
            reads.to_str().unwrap(),
            "-r",
            reference.to_str().unwrap(),
            "-o",
            out_dir.to_str().unwrap(),
            "--prefix",
            prefix,
            "--threads",
            "1",
            "--batch-size",
            "1",
            "--batch-rounds",
            "4",
        ])
        .status()
        .expect("rerun flow overlap mode");
    assert!(status.success());
    assert!(per_gene_unused.exists());

    let summary = fs::read_to_string(out_dir.join("cluster_batch_summary.txt"))
        .expect("read cluster summary");
    assert!(summary.lines().any(|line| line == "processed\t1"));
}
