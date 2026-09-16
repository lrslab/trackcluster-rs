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

fn bed_names(path: &Path) -> Vec<String> {
    fs::read_to_string(path)
        .expect("read BED")
        .lines()
        .filter_map(|line| line.split('\t').nth(3).map(ToOwned::to_owned))
        .collect()
}

fn terminal_bed(name: &str, start: u32, end: u32, score: u32, reverse: bool) -> String {
    let (start, end, sizes, starts, strand) = if reverse {
        (
            2000 - end,
            2000 - start,
            format!("{},{}", end - 500, 200 - start),
            format!("0,{}", end - 200),
            '-',
        )
    } else {
        (
            start,
            end,
            format!("{},{}", 200 - start, end - 500),
            format!("0,{}", 500 - start),
            '+',
        )
    };
    format!("chr1\t{start}\t{end}\t{name}\t{score}\t{strand}\t0\t0\t0\t2\t{sizes},\t{starts},\n")
}

#[test]
fn clusterj_cli_mixed_score_duplicates_are_order_and_batch_invariant() {
    let root = fresh_temp_dir("clusterj_mixed_score_duplicates");
    let reads = root.join("reads.bed");
    let reference = root.join("reference.bed");
    fs::write(&reference, terminal_bed("ref", 0, 1000, 100, false)).unwrap();
    let records = [
        terminal_bed("sl1", 100, 1000, 12, false),
        terminal_bed("sl2", 100, 1000, 12, false),
        terminal_bed("ordinary", 100, 1000, 0, false),
    ];
    let mut expected = None;
    for order in [
        [0, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ] {
        fs::write(&reads, order.map(|index| records[index].as_str()).concat()).unwrap();
        for batch in [0, 1, 500] {
            let out = root.join("isoforms.bed");
            let output = Command::new(env!("CARGO_BIN_EXE_trackcluster"))
                .args(["clusterj", "--reads"])
                .arg(&reads)
                .arg("--reference")
                .arg(&reference)
                .arg("--out")
                .arg(&out)
                .args(["--sw-score", "11", "--batch-size", &batch.to_string()])
                .output()
                .unwrap();
            assert_success(&output, "clusterj mixed-score duplicates");
            let isoforms = trackcluster_rs::io::bed::read_bed12(&out)
                .unwrap()
                .collect::<Result<Vec<_>, _>>()
                .unwrap();
            assert_eq!(isoforms.len(), 2, "order={order:?} batch={batch}");
            let alt = isoforms.iter().find(|tx| tx.tx_start.get() == 100).unwrap();
            assert_eq!(alt.tx_end.get(), 1000);
            assert_eq!(alt.score, 12);
            let mut mapping = normalized_lines(&out.with_extension("read_to_isoform.tsv"));
            mapping.sort();
            assert_eq!(
                mapping,
                ["ordinary", "sl1", "sl2"].map(|name| format!("{name}\t{}", alt.name))
            );
            assert!(normalized_lines(&out.with_extension("unused.bed")).is_empty());
            let mut bed = normalized_lines(&out);
            bed.sort();
            let result = (bed, mapping);
            assert_eq!(expected.get_or_insert_with(|| result.clone()), &result);
        }
    }
}

#[test]
fn clusterj_cli_retains_supported_three_prime_end_with_a_singleton_intermediate() {
    let root = fresh_temp_dir("clusterj_three_prime_intermediate");
    let reads = root.join("reads.bed");
    let reference = root.join("reference.bed");
    fs::write(&reference, terminal_bed("ref", 0, 1000, 100, false)).unwrap();
    let supported: String = (0..5)
        .map(|index| terminal_bed(&format!("alt{index}"), 100, 800, 0, false))
        .collect();
    for bridge in [false, true] {
        let mut input = supported.clone();
        if bridge {
            input.push_str(&terminal_bed("singleton", 100, 830, 0, false));
        }
        fs::write(&reads, input).unwrap();
        for batch in [0, 2, 500] {
            let out = root.join("isoforms.bed");
            let output = Command::new(env!("CARGO_BIN_EXE_trackcluster"))
                .args(["clusterj", "--reads"])
                .arg(&reads)
                .arg("--reference")
                .arg(&reference)
                .arg("--out")
                .arg(&out)
                .args(["--batch-size", &batch.to_string()])
                .output()
                .unwrap();
            assert_success(&output, "clusterj 3prime singleton intermediate");
            let isoforms = trackcluster_rs::io::bed::read_bed12(&out)
                .unwrap()
                .collect::<Result<Vec<_>, _>>()
                .unwrap();
            assert_eq!(isoforms.len(), 2, "bridge={bridge} batch={batch}");
            let alt = isoforms.iter().find(|tx| tx.tx_end.get() == 800).unwrap();
            assert_eq!(alt.tx_start.get(), 100);
            let mapping = normalized_lines(&out.with_extension("read_to_isoform.tsv"));
            for index in 0..5 {
                let pairs: Vec<_> = mapping
                    .iter()
                    .filter(|line| line.starts_with(&format!("alt{index}\t")))
                    .collect();
                assert_eq!(pairs, [&format!("alt{index}\t{}", alt.name)]);
            }
            if bridge {
                assert!(mapping.contains(&"singleton\tref".to_owned()));
            }
            assert!(normalized_lines(&out.with_extension("unused.bed")).is_empty());
        }
    }
}

#[test]
fn clusterj_cli_terminal_evidence_crosses_production_batches_on_both_strands() {
    let root = fresh_temp_dir("clusterj_terminal_production_batches");
    let reads = root.join("reads.bed");
    let reference = root.join("reference.bed");
    let mut refs = String::new();
    let mut input = String::new();
    let mut read_names = std::collections::BTreeSet::new();
    for reverse in [false, true] {
        let label = if reverse { "minus" } else { "plus" };
        refs.push_str(&terminal_bed(
            &format!("ref_{label}"),
            0,
            1000,
            100,
            reverse,
        ));
        // Distinct structures keep this larger than a batch after exact-duplicate coalescing.
        for index in 0..502 {
            let name = format!("full_{label}_{index}");
            input.push_str(&terminal_bed(
                &name,
                index % 100,
                995 + index / 100,
                0,
                reverse,
            ));
            read_names.insert(name);
            if index >= 497 {
                let name = format!("alt_{label}_{index}");
                input.push_str(&terminal_bed(&name, 100, 800, 0, reverse));
                read_names.insert(name);
            }
        }
        let name = format!("singleton_{label}");
        input.push_str(&terminal_bed(&name, 100, 830, 0, reverse));
        read_names.insert(name);
    }
    fs::write(&reference, refs).unwrap();
    fs::write(&reads, input).unwrap();
    let mut expected = None;
    for batch in [499, 500, 501] {
        for threads in [1, 2] {
            let out = root.join("isoforms.bed");
            let output = Command::new(env!("CARGO_BIN_EXE_trackcluster"))
                .args(["clusterj", "--reads"])
                .arg(&reads)
                .arg("--reference")
                .arg(&reference)
                .arg("--out")
                .arg(&out)
                .args([
                    "--batch-size",
                    &batch.to_string(),
                    "--threads",
                    &threads.to_string(),
                ])
                .output()
                .unwrap();
            assert_success(&output, "clusterj production terminal batches");
            let isoforms = trackcluster_rs::io::bed::read_bed12(&out)
                .unwrap()
                .collect::<Result<Vec<_>, _>>()
                .unwrap();
            assert_eq!(isoforms.len(), 4, "batch={batch} threads={threads}");
            let mut mapping = normalized_lines(&out.with_extension("read_to_isoform.tsv"));
            mapping.sort();
            let represented: std::collections::BTreeSet<_> = mapping
                .iter()
                .map(|line| line.split_once('\t').unwrap().0.to_owned())
                .collect();
            assert_eq!(represented, read_names);
            for (label, start, end) in [("plus", 100, 800), ("minus", 1200, 1900)] {
                let alt = isoforms
                    .iter()
                    .find(|tx| (tx.tx_start.get(), tx.tx_end.get()) == (start, end))
                    .unwrap();
                for index in 497..502 {
                    let name = format!("alt_{label}_{index}");
                    let targets: Vec<_> = mapping
                        .iter()
                        .filter(|line| line.starts_with(&format!("{name}\t")))
                        .collect();
                    assert_eq!(targets, [&format!("{name}\t{}", alt.name)]);
                }
            }
            assert!(normalized_lines(&out.with_extension("unused.bed")).is_empty());
            let mut bed = normalized_lines(&out);
            bed.sort();
            let result = (bed, mapping);
            assert_eq!(expected.get_or_insert_with(|| result.clone()), &result);
        }
    }
}

#[test]
fn clusterj_matches_golden_outputs() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");

    let reads = repo_path("tests/fixtures/reads.bed");
    let reference = repo_path("tests/fixtures/ref.bed");
    let golden_isoforms = repo_path("tests/golden/clusterj/isoform.bed");
    let golden_mapping = repo_path("tests/golden/clusterj/isoform.read_to_isoform.tsv");
    let golden_unused = repo_path("tests/golden/clusterj/isoform.unused.bed");

    let out_dir = fresh_temp_dir("clusterj");
    let out_bed = out_dir.join("isoform.bed");

    let output = Command::new(exe)
        .args([
            "clusterj",
            "-s",
            reads.to_str().unwrap(),
            "-r",
            reference.to_str().unwrap(),
            "-o",
            out_bed.to_str().unwrap(),
        ])
        .output()
        .expect("run clusterj");
    assert_success(&output, "clusterj golden run");

    let produced_isoforms = out_bed.clone();
    let produced_mapping = out_bed.with_extension("read_to_isoform.tsv");
    let produced_unused = out_bed.with_extension("unused.bed");

    assert_eq!(
        normalized_lines(&produced_isoforms),
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

#[test]
fn single_gene_cluster_commands_reject_an_output_that_aliases_an_input() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");
    let root = fresh_temp_dir("single_gene_output_alias");
    let reads = root.join("reads.bed");
    let reference = root.join("reference.bed");
    fs::copy(repo_path("tests/fixtures/reads.bed"), &reads).unwrap();
    fs::copy(repo_path("tests/fixtures/ref.bed"), &reference).unwrap();
    let original = fs::read(&reads).unwrap();

    for command in ["clusterj", "cluster"] {
        let output = Command::new(exe)
            .args([command, "--reads"])
            .arg(&reads)
            .arg("--reference")
            .arg(&reference)
            .arg("--out")
            .arg(&reads)
            .output()
            .unwrap();
        assert!(
            !output.status.success(),
            "{command} accepted an input alias"
        );
        assert!(String::from_utf8_lossy(&output.stderr).contains("refer to the same file"));
        assert_eq!(fs::read(&reads).unwrap(), original);
    }
}

#[test]
fn clusterj_plain_bed_reference_is_protected_and_unmatched_reads_are_auditable() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");
    let reads = repo_path("tests/fixtures/cluster_overlap/plain_reads.bed");
    let reference = repo_path("tests/fixtures/cluster_overlap/plain_ref.bed");
    let out_dir = fresh_temp_dir("clusterj_plain_bed");
    let out_bed = out_dir.join("isoform.bed");

    let output = Command::new(exe)
        .args([
            "clusterj",
            "-s",
            reads.to_str().unwrap(),
            "-r",
            reference.to_str().unwrap(),
            "-o",
            out_bed.to_str().unwrap(),
        ])
        .output()
        .expect("run clusterj with plain BED12 inputs");

    assert!(
        output.status.success(),
        "clusterj failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(bed_names(&out_bed), vec!["ref_plain"]);
    assert_eq!(
        normalized_lines(&out_bed.with_extension("read_to_isoform.tsv")),
        vec!["read_match\tref_plain"]
    );
    assert_eq!(
        bed_names(&out_bed.with_extension("unused.bed")),
        vec!["read_wrong_strand", "read_disjoint", "read_wrong_chrom"]
    );

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("input_reads=4"),
        "missing summary: {stderr}"
    );
    assert!(
        stderr.contains("represented_reads=1"),
        "missing summary: {stderr}"
    );
    assert!(
        stderr.contains("mapping_rows=1"),
        "missing summary: {stderr}"
    );
    assert!(stderr.contains("rare_reads=0"), "missing summary: {stderr}");
    assert!(
        stderr.contains("unmatched_reads=3"),
        "missing summary: {stderr}"
    );
    assert!(
        stderr.contains("unused_reads=3"),
        "missing summary: {stderr}"
    );
}

#[test]
fn clusterj_rejects_correction_that_would_create_an_empty_exon() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");
    let input_dir = fresh_temp_dir("clusterj_invalid_snap_input");
    let reads = input_dir.join("reads.bed");
    let reference = input_dir.join("reference.bed");
    fs::write(
        &reference,
        "chr1\t80\t220\tref\t0\t+\t80\t220\t0\t2\t10,20,\t0,120,\n",
    )
    .expect("write correction reference");
    fs::write(
        &reads,
        "chr1\t100\t210\tread\t100\t+\t100\t210\t0\t2\t1,9,\t0,101,\n",
    )
    .expect("write correction read");
    let out_dir = fresh_temp_dir("clusterj_invalid_snap_output");
    let out_bed = out_dir.join("isoform.bed");

    let output = Command::new(exe)
        .args([
            "clusterj",
            "-s",
            reads.to_str().unwrap(),
            "-r",
            reference.to_str().unwrap(),
            "-o",
            out_bed.to_str().unwrap(),
            "--junction-correction-offset",
            "15",
            "--junction-correction-min-support",
            "5",
            "--sw-score",
            "11",
            "--sl-5prime-min-support",
            "1",
            "--sl-same-junction-5prime-offset",
            "0",
        ])
        .output()
        .expect("run clusterj invalid-snap regression");
    assert_success(&output, "clusterj invalid-snap regression");

    let unused = out_bed.with_extension("unused.bed");
    assert_eq!(bed_names(&out_bed), vec!["ref"]);
    assert_eq!(bed_names(&unused), vec!["read"]);
    assert!(String::from_utf8_lossy(&output.stderr).contains("rare_reads=1"));

    for path in [&out_bed, &unused] {
        let validation = Command::new(exe)
            .args(["validate-bed", "--input", path.to_str().unwrap()])
            .output()
            .expect("strictly validate clusterj regression output");
        assert_success(
            &validation,
            "strict validation of clusterj regression output",
        );
    }
}

#[test]
fn single_gene_cluster_commands_skip_only_bad_read_tracks() {
    let exe = env!("CARGO_BIN_EXE_trackcluster");
    let input_dir = fresh_temp_dir("single_gene_bad_read_input");
    let reads = input_dir.join("dirty.bed");
    let good = fs::read_to_string(repo_path("tests/fixtures/reads.bed")).unwrap();
    fs::write(&reads, format!("not-a-bed-record\n{good}")).unwrap();
    let reference = repo_path("tests/fixtures/ref.bed");

    for command in ["clusterj", "cluster"] {
        let out_dir = fresh_temp_dir(&format!("single_gene_{command}_bad_read"));
        let out_bed = out_dir.join("isoform.bed");
        let output = Command::new(exe)
            .args([
                command,
                "-s",
                reads.to_str().unwrap(),
                "-r",
                reference.to_str().unwrap(),
                "-o",
                out_bed.to_str().unwrap(),
            ])
            .output()
            .expect("run single-gene clustering");
        assert!(
            output.status.success(),
            "{command} stopped on one bad read: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        let mapping = fs::read_to_string(out_bed.with_extension("read_to_isoform.tsv")).unwrap();
        assert!(mapping.contains("read_trunc\t"), "{mapping}");
        let rejected = fs::read_to_string(out_bed.with_extension("rejected_reads.tsv")).unwrap();
        assert_eq!(rejected.lines().count(), 2, "{rejected}");
        assert!(String::from_utf8_lossy(&output.stderr).contains("excluded 1 malformed read"));
    }

    let strict_out = input_dir.join("strict.bed");
    let strict = Command::new(exe)
        .args([
            "clusterj",
            "-s",
            reads.to_str().unwrap(),
            "-r",
            reference.to_str().unwrap(),
            "-o",
            strict_out.to_str().unwrap(),
            "--invalid-read-policy",
            "fail",
        ])
        .output()
        .expect("run strict single-gene clustering");
    assert!(!strict.status.success());
    assert!(!strict_out.exists());
}
