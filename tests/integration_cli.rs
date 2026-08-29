use std::process::Command;

#[test]
fn help_and_version_snapshots_are_stable() {
    let trackcluster = env!("CARGO_BIN_EXE_trackcluster");
    let help = Command::new(trackcluster).arg("--help").output().unwrap();
    assert!(help.status.success());
    assert_eq!(
        String::from_utf8(help.stdout).unwrap(),
        include_str!("golden/cli/trackcluster-help.txt")
    );

    let version = Command::new(trackcluster)
        .arg("--version")
        .output()
        .unwrap();
    assert!(version.status.success());
    assert_eq!(
        String::from_utf8(version.stdout).unwrap(),
        include_str!("golden/cli/trackcluster-version.txt")
    );

    let batch = env!("CARGO_BIN_EXE_clusterj_batch");
    let version = Command::new(batch).arg("--version").output().unwrap();
    assert!(version.status.success());
    assert_eq!(
        String::from_utf8(version.stdout).unwrap(),
        include_str!("golden/cli/clusterj-batch-version.txt")
    );
}

#[test]
fn clusterj_batch_help_has_versioned_command_metadata() {
    let output = Command::new(env!("CARGO_BIN_EXE_clusterj_batch"))
        .arg("--help")
        .output()
        .unwrap();
    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(stdout.starts_with("Batch runner for per-gene TrackCluster junction clustering\n"));
    assert!(stdout.contains("Usage: clusterj_batch"));
    assert!(stdout.contains("-V, --version"));
}

#[test]
fn high_expression_gene_cap_defaults_are_consistent() {
    let expected_default = format!(
        "[default: {}]",
        trackcluster_rs::flow::config::DEFAULT_MAX_READS_PER_GENE
    );
    let flow = Command::new(env!("CARGO_BIN_EXE_trackcluster"))
        .args(["flow", "--help"])
        .output()
        .unwrap();
    assert!(flow.status.success());
    let flow_help = String::from_utf8(flow.stdout).unwrap();
    assert!(
        flow_help.contains("--max-reads-per-gene <MAX_READS_PER_GENE>"),
        "{flow_help}"
    );
    assert!(flow_help.contains(&expected_default), "{flow_help}");

    let batch = Command::new(env!("CARGO_BIN_EXE_clusterj_batch"))
        .arg("--help")
        .output()
        .unwrap();
    assert!(batch.status.success());
    let batch_help = String::from_utf8(batch.stdout).unwrap();
    assert!(
        batch_help.contains("--max-reads-per-gene <MAX_READS_PER_GENE>"),
        "{batch_help}"
    );
    assert!(batch_help.contains(&expected_default), "{batch_help}");

    let clusterj = Command::new(env!("CARGO_BIN_EXE_trackcluster"))
        .args(["clusterj", "--help"])
        .output()
        .unwrap();
    assert!(clusterj.status.success());
    let clusterj_help = String::from_utf8(clusterj.stdout).unwrap();
    assert!(
        clusterj_help.contains("--max-reads-per-locus <MAX_READS_PER_LOCUS>"),
        "{clusterj_help}"
    );
    assert!(clusterj_help.contains(&expected_default), "{clusterj_help}");
    assert!(
        clusterj_help.contains("--heartbeat-seconds"),
        "{clusterj_help}"
    );
}
