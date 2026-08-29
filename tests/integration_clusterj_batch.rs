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

fn mapped_gene_key(path: &Path, biological_gene: &str) -> String {
    fs::read_to_string(path)
        .expect("read gene path mapping")
        .lines()
        .filter(|line| !line.starts_with('#') && *line != "gene_id\tpath_key")
        .find_map(|line| {
            let (gene, key) = line.split_once('\t')?;
            (gene == biological_gene).then(|| key.to_owned())
        })
        .expect("biological gene is present in path mapping")
}

fn write_gene_inputs(root: &Path, gene: &str, reads: &str, reference: &str) {
    let gene_dir = root.join(gene);
    fs::create_dir_all(&gene_dir).expect("create gene dir");
    fs::write(gene_dir.join(format!("{gene}_nano.bed")), reads).expect("write gene reads");
    fs::write(gene_dir.join(format!("{gene}_gff.bed")), reference).expect("write gene reference");
}

fn write_mixed_outcome_inputs(root: &Path) {
    let reads = fs::read_to_string(repo_path("tests/fixtures/reads.bed")).expect("read reads");
    let reference =
        fs::read_to_string(repo_path("tests/fixtures/ref.bed")).expect("read reference");

    let missing_dir = root.join("MISSING");
    fs::create_dir_all(&missing_dir).expect("create missing gene dir");
    fs::write(missing_dir.join("MISSING_nano.bed"), &reads).expect("write missing gene reads");
    fs::write(missing_dir.join("run.json"), "{corrupt legacy manifest\n")
        .expect("write corrupt failed-gene manifest");

    write_gene_inputs(root, "EMPTY", "", &reference);

    // Legacy filename-only outputs are deliberately stale without run.json.
    write_gene_inputs(root, "LEGACY", &reads, &reference);
    fs::write(root.join("LEGACY/LEGACY_simple_coveragej.bed"), "")
        .expect("write completed isoforms");
    fs::write(root.join("LEGACY/LEGACY_unused.bed"), "").expect("write completed unused");
    fs::write(root.join("LEGACY/LEGACY_read_to_isoform.tsv"), "").expect("write completed mapping");

    write_gene_inputs(root, "FAILED", &reads, "not-a-reference-record\n");
}

fn run_batch(root: &Path, extra_args: &[&str]) -> std::process::Output {
    let exe = env!("CARGO_BIN_EXE_clusterj_batch");
    let mut command = Command::new(exe);
    command.args([
        "--input-root",
        root.to_str().unwrap(),
        "--output-root",
        root.to_str().unwrap(),
        "--threads",
        "1",
        "--heartbeat-seconds",
        "0",
    ]);
    if !extra_args.contains(&"--gene-list") && !extra_args.contains(&"--prepare-reads") {
        let mut genes = fs::read_dir(root)
            .expect("read test input root")
            .filter_map(|entry| {
                let entry = entry.ok()?;
                entry.file_type().ok()?.is_dir().then(|| entry.file_name())
            })
            .filter_map(|name| name.to_str().map(ToOwned::to_owned))
            .filter(|name| !name.starts_with('.'))
            .collect::<Vec<_>>();
        genes.sort();
        let path = root.join(".test_gene_list.txt");
        let contents = if genes.is_empty() {
            String::new()
        } else {
            format!("{}\n", genes.join("\n"))
        };
        fs::write(&path, contents).expect("write explicit test gene list");
        command.args(["--gene-list", path.to_str().unwrap()]);
    }
    command.args(extra_args);
    command.output().expect("run clusterj_batch")
}

fn assert_resume_decision(root: &Path, action: &str, reason: &str) {
    let summary =
        fs::read_to_string(root.join("clusterj_batch_summary.txt")).expect("read batch summary");
    let expected = format!("resume_decision\tGENEA\t{action}\t{reason}");
    assert!(summary.lines().any(|line| line == expected), "{summary}");
}

#[test]
fn clusterj_batch_keeps_compat_behavior() {
    let exe = env!("CARGO_BIN_EXE_clusterj_batch");

    let reads = repo_path("tests/fixtures/reads.bed");
    let reference = repo_path("tests/fixtures/ref.bed");
    let golden_isoforms = repo_path("tests/golden/clusterj/isoform.bed");
    let golden_unused = repo_path("tests/golden/clusterj/isoform.unused.bed");

    let out_dir = fresh_temp_dir("clusterj_batch");
    let prefix = "sample";

    let output = Command::new(exe)
        .args([
            "--prepare-reads",
            reads.to_str().unwrap(),
            "--prepare-reference",
            reference.to_str().unwrap(),
            "--prepare-prefix",
            prefix,
            "--input-root",
            out_dir.to_str().unwrap(),
            "--output-root",
            out_dir.to_str().unwrap(),
            "--threads",
            "1",
            "--force",
        ])
        .output()
        .expect("run clusterj_batch");
    assert_success(&output, "clusterj batch compatibility run");

    let gene_isoform = out_dir.join("GENEA/GENEA_simple_coveragej.bed");
    let gene_unused = out_dir.join("GENEA/GENEA_unused.bed");
    assert!(gene_isoform.exists());
    assert!(gene_unused.exists());
    assert_eq!(
        normalized_lines(&gene_isoform),
        normalized_lines(&golden_isoforms)
    );
    assert_eq!(
        normalized_lines(&gene_unused),
        normalized_lines(&golden_unused)
    );

    let summary = out_dir.join("clusterj_batch_summary.txt");
    assert!(summary.exists());
    let summary_text = fs::read_to_string(summary).unwrap();
    assert!(summary_text.contains("platform_preset\tgeneric"));
    assert!(summary_text.contains("junction_correction_offset\t10"));
    assert!(summary_text.contains("junction_correction_min_support\t5"));
    assert!(summary_text.contains("same_junction_3prime_offset\t50"));
    assert!(summary_text.contains("3prime_cluster_offset\t10"));
    assert!(summary_text.contains("3prime_min_support\t5"));
    assert!(summary_text.contains("total_genes\t1"));
    assert!(summary_text.contains("errors\t0"));
}

#[test]
fn clusterj_batch_requires_explicit_gene_selection_without_inline_prepare() {
    let root = fresh_temp_dir("clusterj_batch_explicit_gene_selection");
    let reads = fs::read_to_string(repo_path("tests/fixtures/reads.bed")).unwrap();
    let reference = fs::read_to_string(repo_path("tests/fixtures/ref.bed")).unwrap();
    write_gene_inputs(&root, "GENEA", &reads, &reference);

    let output = Command::new(env!("CARGO_BIN_EXE_clusterj_batch"))
        .args([
            "--input-root",
            root.to_str().unwrap(),
            "--output-root",
            root.to_str().unwrap(),
            "--heartbeat-seconds",
            "0",
            "--max-reads-per-gene",
            "0",
        ])
        .output()
        .expect("run batch without an authoritative gene list");
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("--gene-list is required"), "{stderr}");
    assert!(stderr.contains("stale-gene reuse"), "{stderr}");
    assert!(!root.join("GENEA/GENEA_simple_coveragej.bed").exists());
}

#[test]
fn clusterj_batch_rejects_gene_keys_reserved_for_batch_artifacts() {
    let root = fresh_temp_dir("clusterj_batch_reserved_gene_key");
    let reads = fs::read_to_string(repo_path("tests/fixtures/reads.bed")).unwrap();
    let reference = fs::read_to_string(repo_path("tests/fixtures/ref.bed")).unwrap();
    let gene = "clusterj_batch_gene_paths.tsv";
    write_gene_inputs(&root, gene, &reads, &reference.replace("GENEA", gene));

    let output = run_batch(&root, &["--max-reads-per-gene", "0"]);
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("reserved top-level pipeline artifact"),
        "{stderr}"
    );
    assert!(!root
        .join(format!("{gene}/{gene}_simple_coveragej.bed"))
        .exists());
}

#[test]
fn clusterj_batch_inline_prepare_ignores_stale_gene_directories() {
    let root = fresh_temp_dir("clusterj_batch_inline_prepare_stale_gene");
    let reads_path = repo_path("tests/fixtures/reads.bed");
    let reference_path = repo_path("tests/fixtures/ref.bed");
    let reads = fs::read_to_string(&reads_path).unwrap();
    let reference = fs::read_to_string(&reference_path).unwrap();
    write_gene_inputs(&root, "STALE", &reads, &reference);

    let output = Command::new(env!("CARGO_BIN_EXE_clusterj_batch"))
        .args([
            "--prepare-reads",
            reads_path.to_str().unwrap(),
            "--prepare-reference",
            reference_path.to_str().unwrap(),
            "--prepare-prefix",
            "sample",
            "--input-root",
            root.to_str().unwrap(),
            "--output-root",
            root.to_str().unwrap(),
            "--threads",
            "1",
            "--heartbeat-seconds",
            "0",
            "--max-reads-per-gene",
            "0",
            "--force",
        ])
        .output()
        .expect("run inline prepare over a root with stale gene inputs");
    assert_success(&output, "inline prepare stale-gene selection");

    let summary = fs::read_to_string(root.join("clusterj_batch_summary.txt")).unwrap();
    assert!(
        summary.lines().any(|line| line == "total_genes\t1"),
        "{summary}"
    );
    assert!(root.join("GENEA/GENEA_simple_coveragej.bed").exists());
    assert!(!root.join("STALE/STALE_simple_coveragej.bed").exists());
}

#[test]
fn clusterj_batch_manifest_verifies_fingerprint_and_outputs_before_resume() {
    let root = fresh_temp_dir("clusterj_batch_artifact_manifest");
    let reads = fs::read_to_string(repo_path("tests/fixtures/reads.bed")).expect("read reads");
    let reference =
        fs::read_to_string(repo_path("tests/fixtures/ref.bed")).expect("read reference");
    write_gene_inputs(&root, "GENEA", &reads, &reference);

    let first = run_batch(&root, &["--force", "--max-reads-per-gene", "0"]);
    assert!(
        first.status.success(),
        "first batch failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert_resume_decision(&root, "rebuild", "forced");

    let summary = fs::read_to_string(root.join("clusterj_batch_summary.txt"))
        .expect("read source identity from batch summary");
    assert!(summary.lines().any(|line| {
        line == "source_fingerprint\tclean"
            || line
                .strip_prefix("source_fingerprint\tsha256:")
                .is_some_and(|hash| {
                    hash.len() == 64 && hash.bytes().all(|byte| byte.is_ascii_hexdigit())
                })
    }));

    let manifest_path = root.join("GENEA/run.json");
    let manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(&manifest_path).expect("read completion manifest"))
            .expect("parse completion manifest");
    assert_eq!(manifest["schema_version"], 3);
    assert_eq!(manifest["status"], "complete");
    assert_eq!(manifest["request"]["gene"], "GENEA");
    assert_eq!(manifest["request"]["options"]["assignment_mode"], "unique");
    assert_eq!(
        manifest["request"]["options"]["unique_assignment_junction_offset"],
        15
    );
    assert!(manifest["request"]["tool"]["package_version"]
        .as_str()
        .is_some_and(|value| !value.is_empty()));
    assert!(manifest["request"]["tool"]["git_commit"]
        .as_str()
        .is_some_and(|value| !value.is_empty()));
    assert!(manifest["request"]["tool"]["source_fingerprint"]
        .as_str()
        .is_some_and(
            |value| value == "clean" || (value.starts_with("sha256:") && value.len() == 71)
        ));
    assert_eq!(
        manifest["request"]["options"]["invalid_read_policy"],
        "skip"
    );
    assert_eq!(manifest["outputs"].as_array().unwrap().len(), 5);
    for output in manifest["outputs"].as_array().unwrap() {
        assert!(output["sha256"]
            .as_str()
            .is_some_and(|hash| hash.starts_with("sha256:") && hash.len() == 71));
        assert!(output["records"].is_u64());
    }

    let resumed = run_batch(&root, &["--max-reads-per-gene", "0"]);
    assert!(resumed.status.success());
    assert_resume_decision(&root, "reuse", "exact_manifest_match");

    let isoforms = root.join("GENEA/GENEA_simple_coveragej.bed");
    fs::write(&isoforms, "interrupted partial output\n").expect("simulate partial replacement");
    let repaired = run_batch(&root, &["--max-reads-per-gene", "0"]);
    assert!(
        repaired.status.success(),
        "repair failed: {}",
        String::from_utf8_lossy(&repaired.stderr)
    );
    assert_resume_decision(&root, "rebuild", "output_hash_mismatch:isoforms");
    assert!(!fs::read_to_string(&isoforms)
        .unwrap()
        .contains("interrupted partial output"));

    fs::remove_file(&manifest_path).expect("simulate interrupted run before manifest publish");
    fs::write(&isoforms, "interrupted partial output\n").expect("simulate partial output");
    let recovered = run_batch(&root, &["--max-reads-per-gene", "0"]);
    assert!(recovered.status.success());
    assert_resume_decision(&root, "rebuild", "manifest_missing");

    fs::write(root.join("GENEA/GENEA_nano.bed"), format!("{reads}\n"))
        .expect("change reads bytes without changing records");
    let reads_changed = run_batch(&root, &["--max-reads-per-gene", "0"]);
    assert!(reads_changed.status.success());
    assert_resume_decision(&root, "rebuild", "reads_changed");

    fs::write(root.join("GENEA/GENEA_gff.bed"), format!("{reference}\n"))
        .expect("change reference bytes without changing records");
    let reference_changed = run_batch(&root, &["--max-reads-per-gene", "0"]);
    assert!(reference_changed.status.success());
    assert_resume_decision(&root, "rebuild", "reference_changed");

    let enable_downsampling = run_batch(
        &root,
        &["--max-reads-per-gene", "1", "--downsample-seed", "7"],
    );
    assert!(enable_downsampling.status.success());
    assert_resume_decision(&root, "rebuild", "seed_changed");

    let new_seed = run_batch(
        &root,
        &["--max-reads-per-gene", "1", "--downsample-seed", "8"],
    );
    assert!(new_seed.status.success());
    assert_resume_decision(&root, "rebuild", "seed_changed");
}

#[test]
fn clusterj_batch_skips_only_bad_read_tracks_within_a_gene() {
    let root = fresh_temp_dir("clusterj_batch_bad_read_tracks");
    let good = fs::read_to_string(repo_path("tests/fixtures/reads.bed")).expect("read reads");
    let reference =
        fs::read_to_string(repo_path("tests/fixtures/ref.bed")).expect("read reference");
    let dirty = format!(
        "not-a-bed-record\n{good}chr1\t120\t150\tbad_score\tNaN\t+\t120\t150\t0\t1\t30,\t0,\n"
    );
    write_gene_inputs(&root, "GENEA", &dirty, &reference);

    let output = run_batch(&root, &["--force", "--max-reads-per-gene", "0"]);
    assert!(
        output.status.success(),
        "bad read tracks stopped the gene: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let mapping = fs::read_to_string(root.join("GENEA/GENEA_read_to_isoform.tsv"))
        .expect("read retained mapping");
    assert!(mapping.contains("read_trunc\t"), "{mapping}");
    let unused =
        fs::read_to_string(root.join("GENEA/GENEA_unused.bed")).expect("read unused tracks");
    assert!(!unused.contains("bad_score"), "{unused}");
    assert!(!unused.contains("not-a-bed-record"), "{unused}");

    let rejected =
        fs::read_to_string(root.join("GENEA/rejected_reads.tsv")).expect("read rejection report");
    assert_eq!(rejected.lines().count(), 3, "{rejected}");
    assert!(rejected.contains("\t1\t\tparse\t"), "{rejected}");
    assert!(rejected.contains("\t3\tbad_score\tparse\t"), "{rejected}");

    let summary =
        fs::read_to_string(root.join("clusterj_batch_summary.txt")).expect("read summary");
    for expected in [
        "status\tcomplete",
        "invalid_read_policy\tskip",
        "errors\t0",
        "rejected_read_tracks\t2",
        "genes_with_rejected_reads\t1",
        "mergeable_genes\t1",
    ] {
        assert!(summary.lines().any(|line| line == expected), "{summary}");
    }

    let manifest: serde_json::Value = serde_json::from_slice(
        &fs::read(root.join("GENEA/run.json")).expect("read completion manifest"),
    )
    .expect("parse completion manifest");
    let rejected_output = manifest["outputs"]
        .as_array()
        .unwrap()
        .iter()
        .find(|output| output["role"] == "rejected_reads")
        .expect("manifest records rejected-read artifact");
    assert_eq!(rejected_output["records"], 2);

    let count = Command::new(env!("CARGO_BIN_EXE_trackcluster"))
        .args([
            "count",
            "--output-root",
            root.to_str().unwrap(),
            "--prefix",
            "recounted",
            "--reference",
            repo_path("tests/fixtures/ref.bed").to_str().unwrap(),
        ])
        .output()
        .expect("run count-only unique assignment on dirty per-gene reads");
    assert!(
        count.status.success(),
        "unique assignment re-failed on rejected reads: {}",
        String::from_utf8_lossy(&count.stderr)
    );
    assert!(root.join("recounted_isoform_count.csv").exists());

    fs::write(root.join("GENEA/rejected_reads.tsv"), "tampered\n")
        .expect("tamper rejected-read diagnostic");
    let rejected_cache = Command::new(env!("CARGO_BIN_EXE_trackcluster"))
        .args([
            "count",
            "--output-root",
            root.to_str().unwrap(),
            "--prefix",
            "tampered",
            "--reference",
            repo_path("tests/fixtures/ref.bed").to_str().unwrap(),
        ])
        .output()
        .expect("run count-only after diagnostic tamper");
    assert!(!rejected_cache.status.success());
    assert!(!root.join("tampered_isoform.bed").exists());
}

#[test]
fn clusterj_batch_all_bad_reads_publish_verified_empty_result_and_resume() {
    let root = fresh_temp_dir("clusterj_batch_all_bad_read_tracks");
    let reference =
        fs::read_to_string(repo_path("tests/fixtures/ref.bed")).expect("read reference");
    write_gene_inputs(
        &root,
        "GENEA",
        "not-a-bed-record\nchr1\t120\t150\tbad\tNaN\t+\t120\t150\t0\t1\t30,\t0,\n",
        &reference,
    );

    let first = run_batch(&root, &["--max-reads-per-gene", "0"]);
    assert!(
        first.status.success(),
        "all-rejected gene failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    for path in [
        "GENEA/GENEA_simple_coveragej.bed",
        "GENEA/GENEA_unused.bed",
        "GENEA/GENEA_read_to_isoform.tsv",
    ] {
        assert_eq!(fs::metadata(root.join(path)).unwrap().len(), 0, "{path}");
    }
    assert_eq!(
        fs::read_to_string(root.join("GENEA/downsample.tsv"))
            .unwrap()
            .lines()
            .count(),
        1
    );
    assert_eq!(
        fs::read_to_string(root.join("GENEA/rejected_reads.tsv"))
            .unwrap()
            .lines()
            .count(),
        3
    );
    let summary = fs::read_to_string(root.join("clusterj_batch_summary.txt")).unwrap();
    assert!(summary.lines().any(|line| line == "status\tcomplete"));
    assert!(summary
        .lines()
        .any(|line| line == "skipped_no_usable_reads\t1"));
    assert!(summary.lines().any(|line| line == "errors\t0"));

    let diagnostic_before = fs::read(root.join("GENEA/rejected_reads.tsv")).unwrap();
    let manifest_before = fs::read(root.join("GENEA/run.json")).unwrap();
    let resumed = run_batch(&root, &["--max-reads-per-gene", "0"]);
    assert!(resumed.status.success());
    assert_resume_decision(&root, "reuse", "exact_manifest_match");
    let resumed_summary = fs::read_to_string(root.join("clusterj_batch_summary.txt")).unwrap();
    assert!(resumed_summary
        .lines()
        .any(|line| line == "skipped_no_usable_reads\t1"));
    assert!(resumed_summary
        .lines()
        .any(|line| line == "all_reads_rejected_genes\t1"));
    assert_eq!(
        fs::read(root.join("GENEA/rejected_reads.tsv")).unwrap(),
        diagnostic_before
    );
    assert_eq!(
        fs::read(root.join("GENEA/run.json")).unwrap(),
        manifest_before
    );
}

#[test]
fn clusterj_batch_fail_policy_preserves_strict_read_behavior() {
    let root = fresh_temp_dir("clusterj_batch_strict_bad_read");
    let good = fs::read_to_string(repo_path("tests/fixtures/reads.bed")).expect("read reads");
    let reference =
        fs::read_to_string(repo_path("tests/fixtures/ref.bed")).expect("read reference");
    write_gene_inputs(
        &root,
        "GENEA",
        &format!("not-a-bed-record\n{good}"),
        &reference,
    );

    let output = run_batch(
        &root,
        &["--invalid-read-policy", "fail", "--max-reads-per-gene", "0"],
    );
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("no verified gene result"), "{stderr}");
    let summary = fs::read_to_string(root.join("clusterj_batch_summary.txt")).unwrap();
    assert!(summary
        .lines()
        .any(|line| line == "invalid_read_policy\tfail"));
    assert!(summary.lines().any(|line| line == "errors\t1"));
}

#[test]
fn clusterj_batch_downsampling_ignores_rejected_tracks_in_count_and_rng() {
    let clean_root = fresh_temp_dir("clusterj_batch_downsample_clean_reads");
    let dirty_root = fresh_temp_dir("clusterj_batch_downsample_dirty_reads");
    let template = fs::read_to_string(repo_path("tests/fixtures/reads.bed")).expect("read reads");
    let reference =
        fs::read_to_string(repo_path("tests/fixtures/ref.bed")).expect("read reference");
    let valid = (1..=4)
        .map(|index| template.replace("read_trunc", &format!("read_{index}")))
        .collect::<String>();
    write_gene_inputs(&clean_root, "GENEA", &valid, &reference);
    write_gene_inputs(
        &dirty_root,
        "GENEA",
        &format!("not-a-bed-record\n{valid}also-not-a-bed-record\n"),
        &reference,
    );

    for root in [&clean_root, &dirty_root] {
        let output = run_batch(
            root,
            &[
                "--force",
                "--max-reads-per-gene",
                "2",
                "--downsample-seed",
                "7",
            ],
        );
        assert!(
            output.status.success(),
            "downsample batch failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(
            stderr
                .contains("clusterj-batch: subsample gene=GENEA original_reads=4 sampled_reads=2"),
            "{stderr}"
        );
    }

    for suffix in [
        "GENEA_simple_coveragej.bed",
        "GENEA_unused.bed",
        "GENEA_read_to_isoform.tsv",
    ] {
        assert_eq!(
            fs::read(clean_root.join("GENEA").join(suffix)).unwrap(),
            fs::read(dirty_root.join("GENEA").join(suffix)).unwrap(),
            "rejected rows changed deterministic sample for {suffix}"
        );
    }
    let downsample = fs::read_to_string(dirty_root.join("GENEA/downsample.tsv")).unwrap();
    assert!(downsample.contains("GENEA\t4\t2\t2\t"), "{downsample}");
    assert_eq!(
        fs::read_to_string(dirty_root.join("GENEA/rejected_reads.tsv"))
            .unwrap()
            .lines()
            .count(),
        3
    );
}

#[test]
fn clusterj_batch_default_subsamples_a_cox1_sized_locus_before_clustering() {
    let root = fresh_temp_dir("clusterj_batch_default_cox1_cap");
    let template = fs::read_to_string(repo_path("tests/fixtures/reads.bed")).expect("read reads");
    let reference =
        fs::read_to_string(repo_path("tests/fixtures/ref.bed")).expect("read reference");
    let cap = trackcluster_rs::flow::config::DEFAULT_MAX_READS_PER_GENE;
    let original_reads = cap + 1;
    let mut reads = String::with_capacity(template.len() * original_reads);
    for index in 0..original_reads {
        reads.push_str(&template.replace("read_trunc", &format!("cox1_read_{index}")));
    }
    write_gene_inputs(&root, "cox1", &reads, &reference);

    let output = run_batch(&root, &["--force"]);
    assert!(
        output.status.success(),
        "default high-expression cap failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    let expected_log = format!(
        "clusterj-batch: subsample gene=cox1 original_reads={original_reads} sampled_reads={cap}"
    );
    assert!(stderr.contains(&expected_log), "{stderr}");

    let downsample = fs::read_to_string(root.join("cox1/downsample.tsv")).unwrap();
    let scale = original_reads as f64 / cap as f64;
    let expected_record = format!("cox1\t{original_reads}\t{cap}\t{scale}\t");
    assert!(downsample.contains(&expected_record), "{downsample}");
    let summary = fs::read_to_string(root.join("clusterj_batch_summary.txt")).unwrap();
    let expected_summary = format!("max_reads_per_gene\t{cap}");
    assert!(
        summary.lines().any(|line| line == expected_summary),
        "{summary}"
    );
}

#[test]
fn clusterj_batch_uses_encoded_paths_for_unicode_and_long_gene_ids() {
    let exe = env!("CARGO_BIN_EXE_clusterj_batch");
    let reads = repo_path("tests/fixtures/reads.bed");
    let reference_template =
        fs::read_to_string(repo_path("tests/fixtures/ref.bed")).expect("read reference fixture");

    for (label, biological_gene) in [
        ("unicode", "基因-α.1".to_owned()),
        ("long", "very-long-biological-gene-".repeat(30)),
    ] {
        let fixture_dir = fresh_temp_dir(&format!("batch_encoded_fixture_{label}"));
        let reference = fixture_dir.join("reference.bed");
        fs::write(
            &reference,
            reference_template.replace("GENEA", &biological_gene),
        )
        .expect("write special-gene reference");
        let out_dir = fresh_temp_dir(&format!("batch_encoded_{label}"));

        let output = Command::new(exe)
            .args([
                "--prepare-reads",
                reads.to_str().unwrap(),
                "--prepare-reference",
                reference.to_str().unwrap(),
                "--prepare-prefix",
                "sample",
                "--input-root",
                out_dir.to_str().unwrap(),
                "--output-root",
                out_dir.to_str().unwrap(),
                "--threads",
                "1",
                "--heartbeat-seconds",
                "0",
                "--max-reads-per-gene",
                "0",
                "--force",
            ])
            .output()
            .expect("run encoded-gene batch");
        assert!(
            output.status.success(),
            "encoded-gene batch failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );

        let key = mapped_gene_key(
            &out_dir.join("clusterj_batch_gene_paths.tsv"),
            &biological_gene,
        );
        assert!(key.len() < 200);
        assert_ne!(key, biological_gene);
        let gene_dir = out_dir.join(&key);
        assert!(gene_dir.join(format!("{key}_nano.bed")).exists());
        assert!(gene_dir
            .join(format!("{key}_simple_coveragej.bed"))
            .exists());
        assert_eq!(
            fs::read_to_string(gene_dir.join(".trackcluster_gene_id")).unwrap(),
            format!("{biological_gene}\n")
        );
    }
}

#[test]
fn clusterj_batch_rejects_traversal_in_user_gene_lists() {
    let root = fresh_temp_dir("clusterj_batch_gene_list_traversal");
    let outside_name = format!("{}_escape", root.file_name().unwrap().to_string_lossy());
    let outside = root.parent().unwrap().join(&outside_name);
    let gene_list = root.join("genes.txt");
    fs::write(&gene_list, format!("../{outside_name}\n")).unwrap();

    let output = run_batch(
        &root,
        &[
            "--gene-list",
            gene_list.to_str().unwrap(),
            "--max-reads-per-gene",
            "0",
        ],
    );
    assert!(output.status.code().is_some_and(|code| code != 0));
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("invalid gene-list entry"), "{stderr}");
    assert!(stderr.contains("path separators"), "{stderr}");
    assert!(!outside.exists());
}

#[cfg(unix)]
#[test]
fn clusterj_batch_rejects_gene_directory_symlinks_outside_root() {
    use std::os::unix::fs::symlink;

    let parent = fresh_temp_dir("clusterj_batch_symlink_escape");
    let root = parent.join("root");
    let outside = parent.join("outside");
    fs::create_dir_all(&root).unwrap();
    fs::create_dir_all(&outside).unwrap();
    let reads = fs::read_to_string(repo_path("tests/fixtures/reads.bed")).unwrap();
    let reference = fs::read_to_string(repo_path("tests/fixtures/ref.bed")).unwrap();
    fs::write(outside.join("GENEA_nano.bed"), reads).unwrap();
    fs::write(outside.join("GENEA_gff.bed"), reference).unwrap();
    symlink(&outside, root.join("GENEA")).unwrap();
    let gene_list = root.join("genes.txt");
    fs::write(&gene_list, "GENEA\n").unwrap();

    let output = run_batch(
        &root,
        &[
            "--gene-list",
            gene_list.to_str().unwrap(),
            "--max-reads-per-gene",
            "0",
        ],
    );
    assert!(output.status.code().is_some_and(|code| code != 0));
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("escapes output root"), "{stderr}");
    assert!(!outside.join("GENEA_simple_coveragej.bed").exists());
    assert!(!outside.join("GENEA_unused.bed").exists());
    assert!(!outside.join("GENEA_read_to_isoform.tsv").exists());
}

#[test]
fn clusterj_batch_rebuilds_downsample_state_across_partial_resume_and_forced_disable() {
    let root = fresh_temp_dir("clusterj_batch_downsample_resume");
    let reads = format!(
        "{}{}",
        fs::read_to_string(repo_path("tests/fixtures/S1.reads.bed")).expect("read S1"),
        fs::read_to_string(repo_path("tests/fixtures/S2.reads.bed")).expect("read S2")
    );
    let reference =
        fs::read_to_string(repo_path("tests/fixtures/ref.bed")).expect("read reference");
    for gene in ["GENEA", "GENEB"] {
        write_gene_inputs(&root, gene, &reads, &reference);
    }

    let first = run_batch(
        &root,
        &[
            "--force",
            "--max-reads-per-gene",
            "1",
            "--downsample-seed",
            "7",
        ],
    );
    assert!(
        first.status.success(),
        "first batch failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );

    let aggregate_path = root.join("clusterj_batch_downsample.tsv");
    let first_aggregate = fs::read_to_string(&aggregate_path).expect("read aggregate state");
    assert!(first_aggregate.contains("\tinput_fingerprint\n"));
    assert!(first_aggregate.contains("GENEA\t2\t1\t2\t"));
    assert!(first_aggregate.contains("GENEB\t2\t1\t2\t"));
    assert_eq!(first_aggregate.matches("fnv1a64:").count(), 2);

    fs::remove_file(root.join("GENEB/GENEB_read_to_isoform.tsv"))
        .expect("make GENEB a partial resume");
    let resumed = run_batch(
        &root,
        &["--max-reads-per-gene", "1", "--downsample-seed", "7"],
    );
    assert!(
        resumed.status.success(),
        "resumed batch failed: {}",
        String::from_utf8_lossy(&resumed.stderr)
    );
    assert_eq!(
        fs::read_to_string(&aggregate_path).expect("read rebuilt aggregate state"),
        first_aggregate
    );
    let resumed_summary =
        fs::read_to_string(root.join("clusterj_batch_summary.txt")).expect("read resume summary");
    assert!(resumed_summary.lines().any(|line| line == "processed\t1"));
    assert!(resumed_summary.lines().any(|line| line == "skipped\t1"));
    assert!(resumed_summary
        .lines()
        .any(|line| line == "skipped_completed_outputs\t1"));

    let forced_without_downsampling = run_batch(&root, &["--force", "--max-reads-per-gene", "0"]);
    assert!(
        forced_without_downsampling.status.success(),
        "forced batch failed: {}",
        String::from_utf8_lossy(&forced_without_downsampling.stderr)
    );
    assert!(!aggregate_path.exists());
    for gene in ["GENEA", "GENEB"] {
        let state = fs::read_to_string(root.join(format!("{gene}/downsample.tsv")))
            .expect("read explicit empty per-gene state");
        assert_eq!(state.lines().count(), 1);
        assert!(root.join(format!("{gene}/run.json")).exists());
    }
}

#[test]
fn clusterj_batch_fails_and_reports_distinct_gene_outcomes() {
    let root = fresh_temp_dir("clusterj_batch_failure_outcomes");
    write_mixed_outcome_inputs(&root);

    let output = run_batch(
        &root,
        &["--max-reads-per-gene", "0", "--strict-gene-errors"],
    );
    assert!(output.status.code().is_some_and(|code| code != 0));
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("2 required gene(s) failed"), "{stderr}");
    assert!(stderr.contains("no merge/count/description stages were run"));

    let summary =
        fs::read_to_string(root.join("clusterj_batch_summary.txt")).expect("read failure summary");
    for expected in [
        "status\tfailed",
        "gene_error_policy\tstrict",
        "processed\t1",
        "skipped\t1",
        "skipped_completed_outputs\t0",
        "skipped_empty_reads\t1",
        "errors\t2",
        "failed_missing_inputs\t1",
        "failed_processing\t1",
        "failed_panics\t0",
        "mergeable_genes\t2",
        "excluded_failed_genes\t2",
    ] {
        assert!(summary.lines().any(|line| line == expected), "{summary}");
    }
    assert!(
        summary
            .lines()
            .any(|line| line == "resume_decision\tLEGACY\trebuild\tmanifest_missing"),
        "{summary}"
    );

    let errors =
        fs::read_to_string(root.join("clusterj_batch_errors.txt")).expect("read detailed errors");
    assert!(errors.contains("MISSING\tmissing_inputs\tmissing required per-gene input"));
    assert!(errors.contains("FAILED\tprocessing\tparse reference"));
    assert!(errors.contains("expected at least 12 columns"));
    assert!(!root.join("clusterj_batch_downsample.tsv").exists());
}

#[test]
fn clusterj_batch_logs_gene_errors_and_continues_by_default() {
    let root = fresh_temp_dir("clusterj_batch_continue_outcomes");
    write_mixed_outcome_inputs(&root);

    let output = run_batch(&root, &["--max-reads-per-gene", "0"]);
    assert!(
        output.status.success(),
        "default tolerant batch failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("warning: 2 gene(s) failed and were excluded"),
        "{stderr}"
    );
    assert!(
        stderr.contains("continuing with 2 verified gene(s)"),
        "{stderr}"
    );

    let summary =
        fs::read_to_string(root.join("clusterj_batch_summary.txt")).expect("read partial summary");
    for expected in [
        "status\tpartial",
        "gene_error_policy\tcontinue",
        "processed\t1",
        "skipped\t1",
        "errors\t2",
        "mergeable_genes\t2",
        "excluded_failed_genes\t2",
        "infrastructure_errors\t0",
    ] {
        assert!(summary.lines().any(|line| line == expected), "{summary}");
    }

    let errors =
        fs::read_to_string(root.join("clusterj_batch_errors.txt")).expect("read partial errors");
    assert!(errors.contains("MISSING\tmissing_inputs\tmissing required per-gene input"));
    assert!(errors.contains("FAILED\tprocessing\tparse reference"));
}

#[test]
fn clusterj_batch_still_fails_when_no_verified_gene_exists() {
    let root = fresh_temp_dir("clusterj_batch_all_failed");
    let reads = fs::read_to_string(repo_path("tests/fixtures/reads.bed")).expect("read reads");
    write_gene_inputs(&root, "FAILED", &reads, "not-a-reference-record\n");

    let output = run_batch(&root, &["--max-reads-per-gene", "0"]);
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("no verified gene result"), "{stderr}");
    let summary =
        fs::read_to_string(root.join("clusterj_batch_summary.txt")).expect("read failed summary");
    assert!(summary.lines().any(|line| line == "status\tfailed"));
    assert!(summary.lines().any(|line| line == "mergeable_genes\t0"));
}

#[test]
fn partial_batch_rebuilds_downsample_state_from_verified_genes_only() {
    let root = fresh_temp_dir("clusterj_batch_partial_downsample");
    let one_read =
        fs::read_to_string(repo_path("tests/fixtures/reads.bed")).expect("read fixture read");
    let two_reads = format!(
        "{}{}",
        one_read,
        one_read.replace("read_trunc", "read_trunc_second")
    );
    let reference =
        fs::read_to_string(repo_path("tests/fixtures/ref.bed")).expect("read reference");
    write_gene_inputs(
        &root,
        "GOOD",
        &two_reads.replace("GENEA", "GOOD"),
        &reference.replace("GENEA", "GOOD"),
    );

    let failed = root.join("FAILED");
    fs::create_dir_all(&failed).expect("create failed gene");
    fs::write(failed.join("FAILED_nano.bed"), &one_read).expect("write failed reads");
    fs::write(
        root.join("clusterj_batch_downsample.tsv"),
        "gene\toriginal_reads\tsampled_reads\tscale_factor\tseed\ttarget_reads\tinput_fingerprint\nFAILED\t20\t10\t2\t1\t10\tfnv1a64:0123456789abcdef\n",
    )
    .expect("write stale failed-gene aggregate");

    let output = run_batch(
        &root,
        &["--max-reads-per-gene", "1", "--downsample-seed", "7"],
    );
    assert!(
        output.status.success(),
        "partial downsample batch failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let aggregate = fs::read_to_string(root.join("clusterj_batch_downsample.tsv"))
        .expect("read rebuilt partial aggregate");
    assert!(aggregate.contains("GOOD\t2\t1\t2\t"), "{aggregate}");
    assert!(!aggregate.contains("FAILED\t"), "{aggregate}");
}
