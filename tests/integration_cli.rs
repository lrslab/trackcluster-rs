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
