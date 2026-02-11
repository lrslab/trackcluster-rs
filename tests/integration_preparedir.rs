use std::fs;
use std::path::PathBuf;
use std::time::{SystemTime, UNIX_EPOCH};

use trackcluster_rs::annotate::addgene::AddGeneOpts;

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

#[test]
fn preparedir_creates_gene_folders_and_files() {
    let repo_root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let reads = repo_root.join("tests/fixtures/reads.bed");
    let reference = repo_root.join("tests/fixtures/ref.bed");

    let out_root = fresh_temp_dir("preparedir");
    let prefix = "sample";

    let res = trackcluster_rs::flow::preparedir::prepare_dir_from_paths(
        &reads,
        &reference,
        &out_root,
        prefix,
        AddGeneOpts::default(),
    )
    .expect("preparedir run");

    assert_eq!(res.genes, vec!["GENEA".to_owned()]);

    let gene_list = out_root.join(format!("{prefix}_gene.txt"));
    assert!(gene_list.exists());
    assert_eq!(fs::read_to_string(gene_list).unwrap(), "GENEA\n");

    let dedup = out_root.join(format!("{prefix}_dedup.bed"));
    assert!(dedup.exists());
    let dedup_reads: Vec<trackcluster_rs::model::Transcript> =
        trackcluster_rs::io::bed::read_bed12(&dedup)
            .unwrap()
            .collect::<Result<Vec<_>, trackcluster_rs::io::bed::BedError>>()
            .unwrap();
    assert_eq!(dedup_reads.len(), 1);
    assert_eq!(dedup_reads[0].name, "read_trunc");
    assert_eq!(
        dedup_reads[0].extra_fields.get(5).map(|s| s.as_str()),
        Some("none")
    );

    let novel = out_root.join(format!("{prefix}_novel.bed"));
    assert!(novel.exists());
    let novel_reads: Vec<trackcluster_rs::model::Transcript> =
        trackcluster_rs::io::bed::read_bed12(&novel)
            .unwrap()
            .collect::<Result<Vec<_>, trackcluster_rs::io::bed::BedError>>()
            .unwrap();
    assert_eq!(novel_reads.len(), 0);

    let gene_dir = out_root.join("GENEA");
    let gene_reads_path = gene_dir.join("GENEA_nano.bed");
    let gene_ref_path = gene_dir.join("GENEA_gff.bed");
    assert!(gene_reads_path.exists());
    assert!(gene_ref_path.exists());

    let gene_reads: Vec<trackcluster_rs::model::Transcript> =
        trackcluster_rs::io::bed::read_bed12(&gene_reads_path)
            .unwrap()
            .collect::<Result<Vec<_>, trackcluster_rs::io::bed::BedError>>()
            .unwrap();
    assert_eq!(gene_reads.len(), 1);
    assert_eq!(
        gene_reads[0].extra_fields.get(5).map(|s| s.as_str()),
        Some("GENEA")
    );

    let gene_refs: Vec<trackcluster_rs::model::Transcript> =
        trackcluster_rs::io::bed::read_bed12(&gene_ref_path)
            .unwrap()
            .collect::<Result<Vec<_>, trackcluster_rs::io::bed::BedError>>()
            .unwrap();
    assert_eq!(gene_refs.len(), 2);
    for tx in &gene_refs {
        assert_eq!(tx.extra_fields.get(5).map(|s| s.as_str()), Some("GENEA"));
    }
}
