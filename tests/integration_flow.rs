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

fn count_sum(path: &Path) -> f64 {
    let mut reader = csv::Reader::from_path(path).expect("read count csv");
    assert_eq!(
        reader
            .headers()
            .expect("count CSV header")
            .iter()
            .collect::<Vec<_>>(),
        ["gene", "isoform_id", "count"]
    );
    reader
        .records()
        .map(|record| {
            record.expect("count CSV row")[2]
                .parse::<f64>()
                .expect("parse count")
        })
        .sum()
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

#[test]
fn flow_rejects_gene_keys_reserved_for_later_merged_outputs() {
    let fixture_dir = fresh_temp_dir("flow_reserved_gene_fixture");
    let reads = repo_path("tests/fixtures/reads.bed");
    let reference = fixture_dir.join("reference.bed");
    let reference_text = fs::read_to_string(repo_path("tests/fixtures/ref.bed")).unwrap();
    fs::write(
        &reference,
        reference_text.replace("GENEA", "sample_isoform.bed"),
    )
    .unwrap();
    let out_dir = fresh_temp_dir("flow_reserved_gene_output");

    let output = Command::new(env!("CARGO_BIN_EXE_trackcluster"))
        .args([
            "flow",
            "--reads",
            reads.to_str().unwrap(),
            "--reference",
            reference.to_str().unwrap(),
            "--output-root",
            out_dir.to_str().unwrap(),
            "--prefix",
            "sample",
            "--threads",
            "1",
            "--force",
        ])
        .output()
        .expect("run flow with a reserved gene path key");
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("reserved top-level pipeline artifact"),
        "{stderr}"
    );
    assert!(!out_dir.join("sample_gene.txt").exists());
    assert!(!out_dir.join("sample_dedup.bed").exists());
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
    fs::write(
        out_dir.join(format!("{prefix}_sqanti_structural_category.tsv")),
        "stale classification\n",
    )
    .expect("write retired classification artifact");
    for suffix in [
        ".mod_join_qc.tsv",
        ".mod_site_join_qc.tsv",
        ".isoform_mod_sites.tsv",
        ".isoform_mod_design.tsv",
        ".isoform_mod_contrasts.tsv",
    ] {
        fs::write(out_dir.join(format!("{prefix}{suffix}")), "stale\n")
            .expect("write stale optional modification artifact");
    }

    let output = Command::new(exe)
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
            "--unique-assignment-junction-offset",
            "7",
            "--force",
        ])
        .output()
        .expect("run flow");
    assert_success(&output, "golden flow");

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
    let provenance =
        fs::read_to_string(out_dir.join(format!("{prefix}_unique_assignment.provenance.tsv")))
            .expect("read unique-assignment provenance");
    assert!(provenance.contains("unique_assignment_junction_offset\t7\n"));
    assert_eq!(count_sum(&count_out), 1.0);

    assert!(out_dir.join("GENEA/GENEA_simple_coveragej.bed").exists());
    assert!(out_dir.join("GENEA/GENEA_unused.bed").exists());
    assert!(out_dir.join("GENEA/GENEA_read_to_isoform.tsv").exists());
    let run_manifest: serde_json::Value = serde_json::from_slice(
        &fs::read(out_dir.join("GENEA/run.json")).expect("read per-gene run manifest"),
    )
    .expect("parse per-gene run manifest");
    assert_eq!(
        run_manifest["request"]["options"]["unique_assignment_junction_offset"],
        7
    );
    let summary =
        fs::read_to_string(out_dir.join("clusterj_batch_summary.txt")).expect("read batch summary");
    assert!(summary.contains("invocation_json\t["));
    assert!(summary.contains("executable_version\t0.3.0\n"));
    assert!(summary.contains("effective_threads\t1\n"));
    assert!(summary.contains("input_sha256\tGENEA\treads\tsha256:"));
    assert!(summary.contains("input_sha256\tGENEA\treference\tsha256:"));

    assert!(out_dir.join(format!("{prefix}_desc.txt")).exists());
    assert!(out_dir.join(format!("{prefix}_class4.txt")).exists());
    assert!(out_dir.join(format!("{prefix}_fusion.txt")).exists());
    assert!(out_dir.join(format!("{prefix}_class12.txt")).exists());
    assert!(!out_dir
        .join(format!("{prefix}_sqanti_structural_category.tsv"))
        .exists());
    for suffix in [
        ".mod_join_qc.tsv",
        ".mod_site_join_qc.tsv",
        ".isoform_mod_sites.tsv",
        ".isoform_mod_design.tsv",
        ".isoform_mod_contrasts.tsv",
    ] {
        assert!(
            !out_dir.join(format!("{prefix}{suffix}")).exists(),
            "stale optional modification artifact survived: {prefix}{suffix}"
        );
    }
    let desc = fs::read_to_string(out_dir.join(format!("{prefix}_desc.txt"))).unwrap();
    assert!(desc.starts_with("#schema\ttrackcluster-description-v2\tdesc\n"));

    let fractional = Command::new(exe)
        .args([
            "flow",
            "--count-only",
            "--assignment-mode",
            "fractional",
            "--unique-assignment-junction-offset",
            "7",
            "-r",
            reference.to_str().unwrap(),
            "-o",
            out_dir.to_str().unwrap(),
            "--prefix",
            prefix,
        ])
        .output()
        .expect("rerun flow with fractional counting");
    assert_success(&fractional, "fractional count-only rerun");
    assert!(!unique_mapping_out.exists());
    assert!(!out_dir
        .join(format!("{prefix}_unique_assignment.provenance.tsv"))
        .exists());
}

#[test]
fn manifest_flow_runs_optional_modification_aggregation_after_unique_assignment() {
    use trackcluster_rs::io::mod_calls::{
        write_assay_metadata_to_writer, write_observations_tsv_to_writer,
    };
    use trackcluster_rs::model::Strand;
    use trackcluster_rs::modification::{
        AssayMetadata, ImplicitSkipPolicy, ModObservation, ModObservationKey, ModSiteKey,
        ObservationState,
    };

    let fixture = fresh_temp_dir("flow_mod_fixture");
    let out_dir = fresh_temp_dir("flow_mod_output");
    fs::write(
        out_dir.join("pooled.isoform_mod_contrasts.tsv"),
        "stale contrast\n",
    )
    .unwrap();
    let reads = repo_path("tests/fixtures/reads.bed");
    let reference = repo_path("tests/fixtures/ref.bed");
    let manifest = fixture.join("samples.tsv");
    fs::write(
        &manifest,
        format!("sample\tgroup\treads\nS1\tcontrol\t{}\n", reads.display()),
    )
    .unwrap();

    let observations_path = fixture.join("S1.observations.tsv");
    write_observations_tsv_to_writer(
        fs::File::create(&observations_path).unwrap(),
        &[ModObservation {
            key: ModObservationKey {
                assay_id: "a1".to_owned(),
                sample: "S1".to_owned(),
                read_id: "S1::read_trunc".to_owned(),
                site: ModSiteKey {
                    chrom: "chr1".to_owned(),
                    pos0: 125,
                    strand: Strand::Plus,
                    mod_code: "A+a".to_owned(),
                },
            },
            probability: Some(0.9),
            observation_state: ObservationState::ExplicitProbability,
            context: None,
            source_transcript_id: None,
            source_pos0: None,
        }],
    )
    .unwrap();
    let assay_path = fixture.join("a1.assay.json");
    write_assay_metadata_to_writer(
        fs::File::create(&assay_path).unwrap(),
        &AssayMetadata {
            schema_version: 1,
            assay_id: "a1".to_owned(),
            caller: "test".to_owned(),
            caller_version: "1".to_owned(),
            model_id: "test".to_owned(),
            chemistry: "RNA004".to_owned(),
            candidate_rule: "all-context-A".to_owned(),
            source_emission_threshold: None,
            source_site_filter: "none".to_owned(),
            candidate_observations_complete: true,
            implicit_skip_policy: ImplicitSkipPolicy::NotApplicable,
            coordinate_source: "synthetic_genomic".to_owned(),
            read_id_mapping: "sample_prefixed".to_owned(),
            source_files: Vec::new(),
        },
    )
    .unwrap();
    let mod_manifest = fixture.join("mod_samples.tsv");
    fs::write(
        &mod_manifest,
        format!(
            concat!(
                "sample\tassay_id\tobservations\tassay_metadata\tcoverage_bam\n",
                "S1\ta1\t{}\t{}\tNA\n"
            ),
            observations_path.display(),
            assay_path.display()
        ),
    )
    .unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_trackcluster"))
        .args(["flow", "--manifest"])
        .arg(&manifest)
        .arg("--reference")
        .arg(&reference)
        .arg("--output-root")
        .arg(out_dir.path())
        .args([
            "--prefix",
            "pooled",
            "--threads",
            "1",
            "--max-reads-per-gene",
            "0",
            "--force",
            "--mod-manifest",
        ])
        .arg(&mod_manifest)
        .args(["--mod-analysis-threshold", "a1=0.5"])
        .output()
        .unwrap();
    assert_success(&output, "flow with modification aggregation");

    let join = fs::read_to_string(out_dir.join("pooled.mod_join_qc.tsv")).unwrap();
    assert!(join.contains("\tS1\t1\t1\t1\t1\t1\t1\t1\t0\t"), "{join}");
    let sites = fs::read_to_string(out_dir.join("pooled.isoform_mod_sites.tsv")).unwrap();
    assert!(sites.contains("\tchr1\t125\t+\tA+a\t"), "{sites}");
    assert!(
        sites.contains("\t1\tNA\t1\t1\t1\t0\t0\t1\t0.9\t"),
        "{sites}"
    );
    assert!(out_dir.join("pooled.isoform_mod_design.tsv").exists());
    assert!(!out_dir.join("pooled.isoform_mod_contrasts.tsv").exists());
    assert!(out_dir.join("pooled_read_to_isoform.unique.tsv").exists());
}

#[test]
fn flow_round_trips_unicode_and_long_gene_ids_through_encoded_paths() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");
    let reads = repo_path("tests/fixtures/reads.bed");
    let reference_template =
        fs::read_to_string(repo_path("tests/fixtures/ref.bed")).expect("read reference fixture");

    for (label, biological_gene) in [
        ("unicode", "基因-α.1".to_owned()),
        ("long", "very-long-biological-gene-".repeat(30)),
    ] {
        let fixture_dir = fresh_temp_dir(&format!("flow_encoded_fixture_{label}"));
        let reference = fixture_dir.join("reference.bed");
        fs::write(
            &reference,
            reference_template.replace("GENEA", &biological_gene),
        )
        .expect("write special-gene reference");

        let out_dir = fresh_temp_dir(&format!("flow_encoded_{label}"));
        let prefix = "sample";
        let output = Command::new(exe)
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
            .output()
            .expect("run encoded-gene flow");
        assert!(
            output.status.success(),
            "encoded-gene flow failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );

        let mapping_path = out_dir.join("sample_gene_paths.tsv");
        let path_key = mapped_gene_key(&mapping_path, &biological_gene);
        assert!(!path_key.contains(['/', '\\']));
        assert!(path_key.len() < 200);
        assert_ne!(path_key, biological_gene);

        let gene_dir = out_dir.join(&path_key);
        assert!(gene_dir.is_dir());
        assert_eq!(
            fs::read_to_string(gene_dir.join(".trackcluster_gene_id")).unwrap(),
            format!("{biological_gene}\n")
        );
        for suffix in [
            "_nano.bed",
            "_gff.bed",
            "_simple_coveragej.bed",
            "_unused.bed",
            "_read_to_isoform.tsv",
        ] {
            assert!(gene_dir.join(format!("{path_key}{suffix}")).exists());
        }
        let batch_mapping = out_dir.join("clusterj_batch_gene_paths.tsv");
        assert_eq!(mapped_gene_key(&batch_mapping, &biological_gene), path_key);
        let summary = fs::read_to_string(out_dir.join("clusterj_batch_summary.txt")).unwrap();
        assert!(summary.contains("gene_path_map\t"));

        // Exercise prefix-scoped path-map recovery rather than relying on the biological gene list.
        fs::remove_file(out_dir.join("sample_gene.txt")).unwrap();
        let count_only = Command::new(exe)
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
            .output()
            .expect("run count-only encoded-gene flow");
        assert!(
            count_only.status.success(),
            "count-only discovery failed: {}",
            String::from_utf8_lossy(&count_only.stderr)
        );
        assert!(out_dir.join("sample_isoform_count.csv").exists());
    }
}

#[test]
fn flow_rejects_traversing_prefix_before_creating_outputs() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");
    let reads = repo_path("tests/fixtures/reads.bed");
    let reference = repo_path("tests/fixtures/ref.bed");
    let out_dir = fresh_temp_dir("flow_traversing_prefix");
    let outside_name = format!("{}_escape", out_dir.file_name().unwrap().to_string_lossy());
    let outside = out_dir.parent().unwrap().join(&outside_name);
    let unsafe_prefix = format!("../{outside_name}");

    let output = Command::new(exe)
        .args([
            "flow",
            "-s",
            reads.to_str().unwrap(),
            "-r",
            reference.to_str().unwrap(),
            "-o",
            out_dir.to_str().unwrap(),
            "--prefix",
            &unsafe_prefix,
            "--threads",
            "1",
        ])
        .output()
        .expect("run flow with unsafe prefix");
    assert!(output.status.code().is_some_and(|code| code != 0));
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("output prefix"), "{stderr}");
    assert!(stderr.contains("path separators"), "{stderr}");
    assert!(!outside.exists());
    assert_eq!(fs::read_dir(&out_dir).unwrap().count(), 0);
}

#[test]
fn flow_count_only_reuses_completed_gene_outputs() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");

    let reads = repo_path("tests/fixtures/reads.bed");
    let reference = repo_path("tests/fixtures/ref.bed");
    let golden_count = repo_path("tests/golden/count/isoform_count.csv");

    let out_dir = fresh_temp_dir("flow_count_only");
    let prefix = "sample";

    let output = Command::new(exe)
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
        .output()
        .expect("run initial flow");
    assert_success(&output, "initial flow before count-only");

    let count_out = out_dir.join(format!("{prefix}_isoform_count.csv"));
    let unique_mapping_out = out_dir.join(format!("{prefix}_read_to_isoform.unique.tsv"));
    fs::remove_file(&count_out).expect("remove count output");
    fs::remove_file(&unique_mapping_out).expect("remove unique mapping output");

    let output = Command::new(exe)
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
        .output()
        .expect("run count-only flow");
    assert_success(&output, "count-only flow");

    assert_eq!(
        normalized_lines(&count_out),
        normalized_lines(&golden_count)
    );
    assert!(unique_mapping_out.exists());
    assert!(out_dir.join(format!("{prefix}_desc.txt")).exists());
}

#[test]
fn flow_count_only_recovers_only_the_prefix_scoped_gene_set() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");
    let first_reads = repo_path("tests/fixtures/reads.bed");
    let first_reference = repo_path("tests/fixtures/ref.bed");
    let second_reads = repo_path("tests/independent/legacy_trackcluster_0_1_8/inputs/reads.bed");
    let second_reference =
        repo_path("tests/independent/legacy_trackcluster_0_1_8/inputs/reference.bed");
    let out_dir = fresh_temp_dir("flow_count_only_prefix_scope");
    let prefix = "scoped";

    for (reads, reference) in [
        (&first_reads, &first_reference),
        (&second_reads, &second_reference),
    ] {
        let output = Command::new(exe)
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
            .output()
            .expect("run prefix-scoped flow");
        assert_success(&output, "prefix-scoped flow");
    }

    let isoform_path = out_dir.join(format!("{prefix}_isoform.bed"));
    let expected = fs::read(&isoform_path).expect("read current merged isoforms");
    let expected_text = String::from_utf8(expected.clone()).expect("merged isoforms are UTF-8");
    assert!(expected_text.contains("legacy_known_ref"));
    assert!(!expected_text.contains("ref_a"));
    assert!(
        out_dir.join("GENEA").is_dir(),
        "stale gene directory remains"
    );

    fs::remove_file(out_dir.join(format!("{prefix}_gene.txt"))).expect("remove current gene list");
    let count_only = Command::new(exe)
        .args([
            "flow",
            "--count-only",
            "-r",
            second_reference.to_str().unwrap(),
            "-o",
            out_dir.to_str().unwrap(),
            "--prefix",
            prefix,
        ])
        .output()
        .expect("run prefix-scoped count-only flow");
    assert_success(&count_only, "prefix-scoped count-only flow");
    assert_eq!(
        fs::read(&isoform_path).expect("read count-only merged isoforms"),
        expected
    );
}

#[test]
fn flow_count_only_rejects_an_empty_prefix_scoped_gene_list() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");
    let reads = repo_path("tests/fixtures/reads.bed");
    let reference = repo_path("tests/fixtures/ref.bed");
    let out_dir = fresh_temp_dir("flow_count_only_empty_gene_list");
    let prefix = "empty";

    let initial = Command::new(exe)
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
        .output()
        .expect("run initial flow before empty gene-list check");
    assert_success(&initial, "initial flow before empty gene-list check");

    let isoform_path = out_dir.join(format!("{prefix}_isoform.bed"));
    let expected = fs::read(&isoform_path).expect("read merged isoforms before rejection");
    fs::write(out_dir.join(format!("{prefix}_gene.txt")), "")
        .expect("truncate prefix-scoped gene list");

    let count_only = Command::new(exe)
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
        .output()
        .expect("run count-only with empty gene list");
    assert!(!count_only.status.success());
    assert!(String::from_utf8_lossy(&count_only.stderr).contains("selected no genes"));
    assert_eq!(
        fs::read(&isoform_path).expect("read preserved merged isoforms"),
        expected
    );
}

#[test]
fn flow_count_only_rejects_tampered_gene_output_before_publishing_merges() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");
    let reads = repo_path("tests/fixtures/reads.bed");
    let reference = repo_path("tests/fixtures/ref.bed");
    let out_dir = fresh_temp_dir("flow_count_only_tampered");
    let prefix = "sample";

    let initial = Command::new(exe)
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
        .output()
        .expect("run initial flow");
    assert_success(&initial, "initial flow before tamper test");

    let gene = fs::read_to_string(out_dir.join(format!("{prefix}_gene.txt")))
        .expect("read gene list")
        .lines()
        .find(|line| !line.trim().is_empty())
        .expect("at least one gene")
        .to_owned();
    let key = mapped_gene_key(&out_dir.join(format!("{prefix}_gene_paths.tsv")), &gene);
    let per_gene_isoforms = out_dir
        .join(&key)
        .join(format!("{key}_simple_coveragej.bed"));
    fs::write(&per_gene_isoforms, "tampered per-gene result\n").expect("tamper per-gene isoforms");

    let merged_paths = [
        out_dir.join(format!("{prefix}_isoform.bed")),
        out_dir.join(format!("{prefix}_unused.bed")),
        out_dir.join(format!("{prefix}_read_to_isoform.tsv")),
        out_dir.join(format!("{prefix}_isoform_count.csv")),
    ];
    for path in &merged_paths {
        fs::write(path, "previous complete generation\n").expect("write sentinel output");
    }

    let rerun = Command::new(exe)
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
        .output()
        .expect("run tampered count-only flow");
    assert!(
        !rerun.status.success(),
        "tampered result was unexpectedly reused"
    );
    let stderr = String::from_utf8_lossy(&rerun.stderr);
    assert!(stderr.contains("stale or unverified"), "{stderr}");
    assert!(stderr.contains("recorded output content"), "{stderr}");
    for path in &merged_paths {
        assert_eq!(
            fs::read_to_string(path).unwrap(),
            "previous complete generation\n"
        );
    }
}

#[test]
fn flow_manifest_writes_multi_sample_usage_outputs() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");

    let manifest = repo_path("tests/fixtures/samples.tsv");
    let reference = repo_path("tests/fixtures/ref.bed");

    let out_dir = fresh_temp_dir("flow_manifest");
    let prefix = "pooled";

    let output = Command::new(exe)
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
        .output()
        .expect("run flow manifest mode");
    assert_success(&output, "manifest flow with pooled reads");

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
        let mut reader = csv::Reader::from_path(path).expect("read aggregate count csv");
        assert_eq!(
            reader
                .headers()
                .expect("aggregate header")
                .iter()
                .collect::<Vec<_>>(),
            ["gene", "isoform_id", "count"]
        );
        for record in reader.records() {
            let fields = record.expect("aggregate row");
            assert_eq!(fields.len(), 3);
            let matrix_total = matrix_totals
                .get(&fields[1])
                .copied()
                .expect("matrix row for isoform");
            let aggregate = fields[2].parse::<f64>().expect("parse aggregate count");
            assert!((aggregate - matrix_total).abs() < 1e-9);
        }
    }
}

#[test]
fn flow_manifest_skips_invalid_tracks_before_sample_tagging() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");
    let input_dir = fresh_temp_dir("flow_manifest_invalid_tracks_input");
    let reads = input_dir.join("sample.bed");
    let good = fs::read_to_string(repo_path("tests/fixtures/reads.bed")).unwrap();
    fs::write(&reads, format!("not-a-bed-record\n{good}")).unwrap();
    let manifest = input_dir.join("samples.tsv");
    fs::write(
        &manifest,
        format!("sample\treads\nS1\t{}\n", reads.display()),
    )
    .unwrap();
    let reference = repo_path("tests/fixtures/ref.bed");
    let out_dir = fresh_temp_dir("flow_manifest_invalid_tracks_output");

    let output = Command::new(exe)
        .args([
            "flow",
            "--manifest",
            manifest.to_str().unwrap(),
            "-r",
            reference.to_str().unwrap(),
            "-o",
            out_dir.to_str().unwrap(),
            "--prefix",
            "pooled",
            "--emit-pooled-reads",
            "--threads",
            "1",
            "--heartbeat-seconds",
            "0",
            "--max-reads-per-gene",
            "0",
            "--force",
        ])
        .output()
        .expect("run manifest flow with invalid read");
    assert!(
        output.status.success(),
        "manifest flow failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let pooled = fs::read_to_string(out_dir.join("pooled_pooled_reads.bed")).unwrap();
    assert!(pooled.contains("\tS1::read_trunc\t"), "{pooled}");
    assert!(!pooled.contains("not-a-bed-record"), "{pooled}");
    let rejected = fs::read_to_string(out_dir.join("pooled_rejected_reads.tsv")).unwrap();
    assert_eq!(rejected.lines().count(), 2, "{rejected}");
    let summary = fs::read_to_string(out_dir.join("clusterj_batch_summary.txt")).unwrap();
    assert!(summary
        .lines()
        .any(|line| line == "prepare_rejected_read_tracks\t1"));
    assert!(summary
        .lines()
        .any(|line| line == "rejected_read_tracks\t1"));
}

#[test]
fn flow_manifest_skips_pooled_reads_by_default() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");

    let manifest = repo_path("tests/fixtures/samples.tsv");
    let reference = repo_path("tests/fixtures/ref.bed");

    let out_dir = fresh_temp_dir("flow_manifest_no_pool");
    let prefix = "pooled";

    let output = Command::new(exe)
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
        .output()
        .expect("run flow manifest mode");
    assert_success(&output, "manifest flow without pooled reads");

    let pooled_reads = out_dir.join(format!("{prefix}_pooled_reads.bed"));
    assert!(!pooled_reads.exists());
}

#[test]
fn flow_manifest_resume_preserves_downsampled_counts_and_usage_byte_for_byte() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");

    let manifest = repo_path("tests/fixtures/samples.tsv");
    let reference = repo_path("tests/fixtures/ref.bed");

    let out_dir = fresh_temp_dir("flow_manifest_downsample");
    let prefix = "pooled";

    let output = Command::new(exe)
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
        .output()
        .expect("run flow manifest downsample");
    assert_success(&output, "manifest flow with downsampling");

    let downsample_state = out_dir.join("clusterj_batch_downsample.tsv");
    assert!(downsample_state.exists());
    let summary_text = fs::read_to_string(&downsample_state).unwrap();
    assert!(summary_text
        .lines()
        .any(|line| line.starts_with("GENEA\t2\t1\t2")));
    assert!(summary_text.contains("input_fingerprint"));
    assert!(summary_text.contains("fnv1a64:"));

    assert!(out_dir.join("GENEA/downsample.tsv").exists());

    let main_count = out_dir.join(format!("{prefix}_isoform_count.csv"));
    let multi_count = out_dir.join(format!("{prefix}.isoform_count.csv"));
    assert!((count_sum(&main_count) - 2.0).abs() < 1e-9);
    assert!((count_sum(&multi_count) - 2.0).abs() < 1e-9);

    let stable_outputs = [
        main_count,
        multi_count,
        out_dir.join(format!("{prefix}.isoform_usage.long.tsv")),
        out_dir.join(format!("{prefix}.isoform_counts.matrix.tsv")),
        out_dir.join(format!("{prefix}.isoform_usage.group.tsv")),
        downsample_state,
    ];
    let before_resume: Vec<Vec<u8>> = stable_outputs
        .iter()
        .map(|path| fs::read(path).expect("read output before resume"))
        .collect();

    let resumed = Command::new(exe)
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
            "--max-reads-per-gene",
            "1",
            "--downsample-seed",
            "1",
        ])
        .output()
        .expect("resume flow manifest downsample");
    assert!(
        resumed.status.success(),
        "resume failed: {}",
        String::from_utf8_lossy(&resumed.stderr)
    );
    for (path, expected) in stable_outputs.iter().zip(before_resume) {
        assert_eq!(
            fs::read(path).expect("read resumed output"),
            expected,
            "{path:?}"
        );
    }
    let batch_summary = fs::read_to_string(out_dir.join("clusterj_batch_summary.txt"))
        .expect("read resumed batch summary");
    assert!(batch_summary.lines().any(|line| line == "processed\t0"));
    assert!(batch_summary.lines().any(|line| line == "skipped\t1"));
    assert!(batch_summary
        .lines()
        .any(|line| line == "skipped_completed_outputs\t1"));
    assert!(batch_summary.lines().any(|line| line == "errors\t0"));
}

#[test]
fn flow_count_only_rebuilds_missing_aggregate_downsample_state() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");
    let manifest = repo_path("tests/fixtures/samples.tsv");
    let reference = repo_path("tests/fixtures/ref.bed");
    let out_dir = fresh_temp_dir("flow_count_only_downsample_recovery");
    let prefix = "downsampled";

    let initial = Command::new(exe)
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
            "--assignment-mode",
            "fractional",
        ])
        .output()
        .expect("run initial downsampled manifest flow");
    assert_success(&initial, "initial downsampled manifest flow");

    let count_path = out_dir.join(format!("{prefix}_isoform_count.csv"));
    let multi_count_path = out_dir.join(format!("{prefix}.isoform_count.csv"));
    let expected_count = fs::read(&count_path).expect("read scaled count output");
    let expected_multi_count = fs::read(&multi_count_path).expect("read scaled multi-count output");
    assert!((count_sum(&count_path) - 2.0).abs() < 1e-9);

    let aggregate = out_dir.join("clusterj_batch_downsample.tsv");
    fs::remove_file(&aggregate).expect("remove aggregate downsample state");
    let count_only = Command::new(exe)
        .args([
            "flow",
            "--count-only",
            "--manifest",
            manifest.to_str().unwrap(),
            "-r",
            reference.to_str().unwrap(),
            "-o",
            out_dir.to_str().unwrap(),
            "--prefix",
            prefix,
            "--assignment-mode",
            "fractional",
        ])
        .output()
        .expect("run count-only after removing aggregate downsample state");
    assert_success(
        &count_only,
        "count-only after removing aggregate downsample state",
    );

    assert_eq!(
        fs::read(&count_path).expect("read recovered count"),
        expected_count
    );
    assert_eq!(
        fs::read(&multi_count_path).expect("read recovered multi-count"),
        expected_multi_count
    );
    let rebuilt = fs::read_to_string(&aggregate).expect("read rebuilt aggregate downsample state");
    assert!(rebuilt
        .lines()
        .any(|line| line.starts_with("GENEA\t2\t1\t2")));
}

#[test]
fn flow_rejects_independent_downsampling_of_multi_gene_molecules() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");
    let input_dir = fresh_temp_dir("flow_multi_gene_downsample_input");
    let reads = input_dir.join("reads.bed");
    let reference = input_dir.join("reference.bed");
    fs::write(
        &reference,
        concat!(
            "chr1\t100\t250\tref_a\t100\t+\t0\t0\t0\t2\t50,50,\t0,100,\tnone\tnone\tnone\t-1,-1,\tisoform_anno\tGENEA\tnone\tnone\n",
            "chr1\t100\t250\tref_b\t100\t+\t0\t0\t0\t2\t50,50,\t0,100,\tnone\tnone\tnone\t-1,-1,\tisoform_anno\tGENEB\tnone\tnone\n"
        ),
    )
    .expect("write overlapping references");
    fs::write(
        &reads,
        concat!(
            "chr1\t100\t250\tr1\t1\t+\t0\t0\t0\t2\t50,50,\t0,100,\tnone\tnone\tnone\t-1,-1,\tnanopore_read\tnone\tnone\tnone\n",
            "chr1\t100\t250\tr2\t1\t+\t0\t0\t0\t2\t50,50,\t0,100,\tnone\tnone\tnone\t-1,-1,\tnanopore_read\tnone\tnone\tnone\n"
        ),
    )
    .expect("write multi-gene reads");

    let safe_out = fresh_temp_dir("flow_multi_gene_no_downsample");
    let safe = Command::new(exe)
        .args([
            "flow",
            "-s",
            reads.to_str().unwrap(),
            "-r",
            reference.to_str().unwrap(),
            "-o",
            safe_out.to_str().unwrap(),
            "--prefix",
            "safe",
            "--threads",
            "1",
            "--force",
            "--assignment-mode",
            "fractional",
            "--max-reads-per-gene",
            "0",
        ])
        .output()
        .expect("run multi-gene flow without downsampling");
    assert_success(&safe, "multi-gene flow without downsampling");
    assert!((count_sum(&safe_out.join("safe_isoform_count.csv")) - 2.0).abs() < 1e-9);

    let unsafe_out = fresh_temp_dir("flow_multi_gene_with_downsample");
    let unsafe_run = Command::new(exe)
        .args([
            "flow",
            "-s",
            reads.to_str().unwrap(),
            "-r",
            reference.to_str().unwrap(),
            "-o",
            unsafe_out.to_str().unwrap(),
            "--prefix",
            "unsafe",
            "--threads",
            "1",
            "--force",
            "--assignment-mode",
            "fractional",
            "--max-reads-per-gene",
            "1",
        ])
        .output()
        .expect("run unsafe multi-gene downsampling regression");
    assert!(!unsafe_run.status.success());
    assert!(String::from_utf8_lossy(&unsafe_run.stderr)
        .contains("independent per-gene downsampling is unsafe"));
    assert!(!unsafe_out.join("unsafe_isoform_count.csv").exists());
}

#[test]
fn flow_single_sample_resume_keeps_scaled_total() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");
    let reads_dir = fresh_temp_dir("flow_single_downsample_reads");
    let reads = reads_dir.join("reads.bed");
    fs::write(
        &reads,
        format!(
            "{}{}",
            fs::read_to_string(repo_path("tests/fixtures/S1.reads.bed")).expect("read S1"),
            fs::read_to_string(repo_path("tests/fixtures/S2.reads.bed")).expect("read S2")
        ),
    )
    .expect("write combined reads");
    let reference = repo_path("tests/fixtures/ref.bed");
    let out_dir = fresh_temp_dir("flow_single_downsample_resume");
    let prefix = "sample";

    for force in [true, false] {
        let mut command = Command::new(exe);
        command.args([
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
            "--max-reads-per-gene",
            "1",
            "--downsample-seed",
            "11",
        ]);
        if force {
            command.arg("--force");
        }
        let output = command.output().expect("run single-sample flow");
        assert!(
            output.status.success(),
            "flow failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert!((count_sum(&out_dir.join("sample_isoform_count.csv")) - 2.0).abs() < 1e-9);
    }
}

#[test]
fn flow_malformed_input_exits_nonzero_before_merged_artifacts() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");
    let input_dir = fresh_temp_dir("flow_malformed_input");
    let reads = input_dir.join("malformed.bed");
    fs::write(&reads, "not-a-bed-record\n").expect("write malformed reads");
    let reference = repo_path("tests/fixtures/ref.bed");
    let out_dir = fresh_temp_dir("flow_malformed_output");
    let prefix = "failed";

    let output = Command::new(exe)
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
        .output()
        .expect("run malformed flow");
    assert!(output.status.code().is_some_and(|code| code != 0));
    assert!(!out_dir.join(format!("{prefix}_isoform.bed")).exists());
    assert!(!out_dir.join(format!("{prefix}_isoform_count.csv")).exists());
    assert!(!out_dir.join(format!("{prefix}_desc.txt")).exists());
}

#[test]
fn flow_skips_bad_raw_read_tracks_and_finishes_with_good_tracks() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");
    let input_dir = fresh_temp_dir("flow_mixed_read_tracks_input");
    let reads = input_dir.join("mixed.bed");
    let good = fs::read_to_string(repo_path("tests/fixtures/reads.bed")).expect("read fixture");
    fs::write(
        &reads,
        format!(
            "not-a-bed-record\n{good}chr1\t120\t150\tbad_score\tNaN\t+\t120\t150\t0\t1\t30,\t0,\n"
        ),
    )
    .expect("write mixed reads");
    let reference = repo_path("tests/fixtures/ref.bed");
    let out_dir = fresh_temp_dir("flow_mixed_read_tracks_output");

    let output = Command::new(exe)
        .args([
            "flow",
            "-s",
            reads.to_str().unwrap(),
            "-r",
            reference.to_str().unwrap(),
            "-o",
            out_dir.to_str().unwrap(),
            "--prefix",
            "mixed",
            "--threads",
            "1",
            "--heartbeat-seconds",
            "0",
            "--max-reads-per-gene",
            "0",
            "--force",
        ])
        .output()
        .expect("run mixed-read flow");
    assert!(
        output.status.success(),
        "mixed-read flow failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!((count_sum(&out_dir.join("mixed_isoform_count.csv")) - 1.0).abs() < 1e-9);
    let rejected = fs::read_to_string(out_dir.join("mixed_rejected_reads.tsv"))
        .expect("read prepare rejection report");
    assert_eq!(rejected.lines().count(), 3, "{rejected}");
    let per_gene_rejected = fs::read_to_string(out_dir.join("GENEA/rejected_reads.tsv"))
        .expect("read per-gene rejection report");
    assert_eq!(per_gene_rejected.lines().count(), 1, "{per_gene_rejected}");
    let summary =
        fs::read_to_string(out_dir.join("clusterj_batch_summary.txt")).expect("read flow summary");
    assert!(summary
        .lines()
        .any(|line| line == "prepare_rejected_read_tracks\t2"));
    assert!(summary
        .lines()
        .any(|line| line == "rejected_read_tracks\t2"));
    assert!(summary.lines().any(|line| line == "errors\t0"));

    let strict_out = fresh_temp_dir("flow_mixed_read_tracks_strict_output");
    let strict = Command::new(exe)
        .args([
            "flow",
            "-s",
            reads.to_str().unwrap(),
            "-r",
            reference.to_str().unwrap(),
            "-o",
            strict_out.to_str().unwrap(),
            "--prefix",
            "strict",
            "--invalid-read-policy",
            "fail",
            "--threads",
            "1",
            "--heartbeat-seconds",
            "0",
        ])
        .output()
        .expect("run strict mixed-read flow");
    assert!(!strict.status.success());
    assert!(!strict_out.join("strict_isoform.bed").exists());
}

#[test]
fn flow_excludes_failed_genes_and_finishes_partial_outputs_by_default() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");
    let input_dir = fresh_temp_dir("flow_partial_input");
    let reads = input_dir.join("reads.bed");
    let reference = input_dir.join("reference.bed");
    fs::write(
        &reads,
        concat!(
            "chr1\t110\t190\tread_good\t0\t+\t110\t190\t0\t1\t80,\t0,\n",
            "chr2\t110\t190\tread_bad\t0\t+\t110\t190\t0\t1\t80,\t0,\n",
        ),
    )
    .expect("write partial-flow reads");
    fs::write(
        &reference,
        concat!(
            "chr1\t100\t200\tref_good\t0\t+\t100\t200\t0\t1\t100,\t0,\tnone\tnone\tnone\tnone\tisoform_anno\tGOOD\n",
            "chr2\t100\t200\tdup_bad\t0\t+\t100\t200\t0\t1\t100,\t0,\tnone\tnone\tnone\tnone\tisoform_anno\tBAD\n",
            "chr2\t100\t220\tdup_bad\t0\t+\t100\t220\t0\t1\t120,\t0,\tnone\tnone\tnone\tnone\tisoform_anno\tBAD\n",
        ),
    )
    .expect("write partial-flow reference");

    let out_dir = fresh_temp_dir("flow_partial_output");
    let stale_bad_dir = out_dir.join("BAD");
    fs::create_dir_all(&stale_bad_dir).expect("create stale BAD folder");
    fs::write(
        stale_bad_dir.join("BAD_simple_coveragej.bed"),
        "chr2\t100\t200\tstale_bad\t0\t+\t100\t200\t0\t1\t100,\t0,\tnone\tnone\tnone\tnone\tisoform_anno\tBAD\n",
    )
    .expect("write stale failed-gene isoform");
    fs::write(stale_bad_dir.join("BAD_unused.bed"), "").expect("write stale unused");
    fs::write(
        stale_bad_dir.join("BAD_read_to_isoform.tsv"),
        "read_bad\tstale_bad\n",
    )
    .expect("write stale mapping");
    fs::write(
        out_dir.join("clusterj_batch_downsample.tsv"),
        "gene\toriginal_reads\tsampled_reads\tscale_factor\tseed\ttarget_reads\tinput_fingerprint\nBAD\t10\t5\t2\t1\t5\tfnv1a64:0123456789abcdef\n",
    )
    .expect("write stale failed-gene downsample state");

    let output = Command::new(exe)
        .args([
            "flow",
            "--reads",
            reads.to_str().unwrap(),
            "--reference",
            reference.to_str().unwrap(),
            "--output-root",
            out_dir.to_str().unwrap(),
            "--prefix",
            "partial",
            "--threads",
            "1",
            "--heartbeat-seconds",
            "0",
            "--max-reads-per-gene",
            "0",
            "--force",
        ])
        .output()
        .expect("run tolerant partial flow");
    assert_success(&output, "tolerant partial flow");

    let summary = fs::read_to_string(out_dir.join("clusterj_batch_summary.txt"))
        .expect("read partial flow summary");
    assert!(summary.lines().any(|line| line == "status\tpartial"));
    assert!(summary.lines().any(|line| line == "mergeable_genes\t1"));
    assert!(summary
        .lines()
        .any(|line| line == "excluded_failed_genes\t1"));
    let errors = fs::read_to_string(out_dir.join("clusterj_batch_errors.txt"))
        .expect("read partial flow errors");
    assert!(errors.contains("BAD\tprocessing"), "{errors}");
    assert!(
        errors.contains("duplicate reference isoform id"),
        "{errors}"
    );

    let merged = fs::read_to_string(out_dir.join("partial_isoform.bed"))
        .expect("read partial merged isoforms");
    let merged_ids: Vec<&str> = merged
        .lines()
        .filter_map(|line| line.split('\t').nth(3))
        .collect();
    assert_eq!(merged_ids, ["ref_good"]);
    assert!(!merged.contains("stale_bad"));
    let count = fs::read_to_string(out_dir.join("partial_isoform_count.csv"))
        .expect("read partial count output");
    assert!(count.contains("GOOD,ref_good,1"), "{count}");
    assert!(!count.contains("BAD"), "{count}");
    assert!(out_dir.join("partial_desc.txt").is_file());
    assert!(!out_dir.join("clusterj_batch_downsample.tsv").exists());

    let strict_out = fresh_temp_dir("flow_partial_strict_output");
    let strict = Command::new(exe)
        .args([
            "flow",
            "--reads",
            reads.to_str().unwrap(),
            "--reference",
            reference.to_str().unwrap(),
            "--output-root",
            strict_out.to_str().unwrap(),
            "--prefix",
            "strict",
            "--threads",
            "1",
            "--heartbeat-seconds",
            "0",
            "--max-reads-per-gene",
            "0",
            "--force",
            "--strict-gene-errors",
        ])
        .output()
        .expect("run strict partial flow");
    assert!(!strict.status.success());
    assert!(!strict_out.join("strict_isoform.bed").exists());
    assert!(!strict_out.join("strict_isoform_count.csv").exists());
    assert!(!strict_out.join("strict_desc.txt").exists());
}

#[test]
fn flow_overlap_mode_runs_end_to_end() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");

    let reads_dir = fresh_temp_dir("flow_overlap_reads");
    let read1 = fs::read_to_string(repo_path("tests/fixtures/S1.reads.bed")).expect("read S1");
    let read2 = fs::read_to_string(repo_path("tests/fixtures/S2.reads.bed")).expect("read S2");
    let reads = reads_dir.join("reads.bed");
    fs::write(&reads, format!("{read1}{read2}")).expect("write overlap reads");
    let reference = repo_path("tests/fixtures/ref.bed");

    let out_dir = fresh_temp_dir("flow_overlap");
    let prefix = "sample";

    let output = Command::new(exe)
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
        .output()
        .expect("run flow overlap mode");
    assert_success(&output, "overlap-mode flow");

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
    let output = Command::new(exe)
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
        .output()
        .expect("rerun flow overlap mode");
    assert_success(&output, "resumed overlap-mode flow");
    assert!(per_gene_unused.exists());

    let summary = fs::read_to_string(out_dir.join("cluster_batch_summary.txt"))
        .expect("read cluster summary");
    assert!(summary.lines().any(|line| line == "processed\t1"));
}
