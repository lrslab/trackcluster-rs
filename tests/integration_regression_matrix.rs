mod common;

use std::fs;
use std::path::Path;
use std::process::{Command, Output};

use common::{assert_success, TestDir};

fn reference_bed_line(chrom: &str, start: u32, name: &str, strand: char, gene: &str) -> String {
    let end = start + 30;
    format!(
        "{chrom}\t{start}\t{end}\t{name}\t100\t{strand}\t{start}\t{end}\t0\t2\t10,10,\t0,20,\tnone\tnone\tnone\t-1,-1,\tisoform_anno\t{gene}\tgroup1\tnone\n"
    )
}

fn read_bed_line(chrom: &str, start: u32, name: &str, strand: char) -> String {
    let end = start + 30;
    format!(
        "{chrom}\t{start}\t{end}\t{name}\t1\t{strand}\t{start}\t{end}\t0\t2\t10,10,\t0,20,\tnone\tnone\tnone\t-1,-1,\tnanopore_read\tnone\tgroup1\tnone\n"
    )
}

fn one_exon_bed_line(
    chrom: &str,
    start: u32,
    end: u32,
    name: &str,
    strand: char,
    record_type: &str,
    gene: &str,
) -> String {
    let size = end - start;
    format!(
        "{chrom}\t{start}\t{end}\t{name}\t1\t{strand}\t{start}\t{end}\t0\t1\t{size},\t0,\tnone\tnone\tnone\t-1,\t{record_type}\t{gene}\tgroup1\tnone\n"
    )
}

fn run_trackcluster(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_trackcluster"))
        .args(args)
        .output()
        .expect("run trackcluster")
}

fn nonempty_lines(path: &Path) -> Vec<String> {
    fs::read_to_string(path)
        .unwrap_or_else(|error| panic!("read {path:?}: {error}"))
        .lines()
        .filter(|line| !line.is_empty())
        .map(ToOwned::to_owned)
        .collect()
}

fn bed_names(path: &Path) -> Vec<String> {
    trackcluster_rs::io::bed::read_bed12(path)
        .unwrap_or_else(|error| panic!("open BED {path:?}: {error}"))
        .map(|record| record.expect("parse BED record").name)
        .collect()
}

fn count_sum(path: &Path) -> f64 {
    let mut reader = csv::Reader::from_path(path).expect("open count CSV");
    assert_eq!(
        reader
            .headers()
            .expect("count header")
            .iter()
            .collect::<Vec<_>>(),
        ["gene", "isoform_id", "count"]
    );
    reader
        .records()
        .map(|row| {
            row.expect("count row")[2]
                .parse::<f64>()
                .expect("numeric count")
        })
        .sum()
}

#[test]
fn flow_conserves_reads_across_multiple_genes_chromosomes_and_strands() {
    let fixture = TestDir::new("regression-multilocus");
    let reads = fixture.path().join("reads.bed");
    let references = fixture.path().join("references.bed");
    let output_root = fixture.path().join("out");

    // Deliberately scramble genomic and gene order. Published output order is part of this test.
    fs::write(
        &references,
        [
            reference_bed_line("chr3", 700, "ref_unknown", '.', "GENE_UNKNOWN"),
            reference_bed_line("chr1", 100, "ref_plus", '+', "GENE_PLUS"),
            reference_bed_line("chr2", 400, "ref_minus", '-', "GENE_MINUS"),
        ]
        .concat(),
    )
    .expect("write multi-locus references");
    fs::write(
        &reads,
        [
            read_bed_line("chr3", 700, "read_unknown", '.'),
            // Same locus as ref_plus but a different strand: it must remain auditable as novel.
            read_bed_line("chr1", 100, "read_wrong_strand", '-'),
            read_bed_line("chr2", 400, "read_minus", '-'),
            read_bed_line("chr1", 100, "read_plus", '+'),
        ]
        .concat(),
    )
    .expect("write multi-locus reads");

    let output = run_trackcluster(&[
        "flow",
        "--reads",
        reads.to_str().unwrap(),
        "--reference",
        references.to_str().unwrap(),
        "--output-root",
        output_root.to_str().unwrap(),
        "--prefix",
        "matrix",
        "--threads",
        "2",
        "--max-reads-per-gene",
        "0",
        "--force",
    ]);
    assert_success(&output, "multi-locus flow");

    assert_eq!(
        nonempty_lines(&output_root.join("matrix_gene.txt")),
        ["GENE_MINUS", "GENE_PLUS", "GENE_UNKNOWN"]
    );
    for gene in ["GENE_MINUS", "GENE_PLUS", "GENE_UNKNOWN"] {
        assert!(output_root.join(gene).is_dir(), "missing folder for {gene}");
    }

    assert_eq!(
        bed_names(&output_root.join("matrix_isoform.bed")),
        ["ref_minus", "ref_plus", "ref_unknown"]
    );
    let mappings = nonempty_lines(&output_root.join("matrix_read_to_isoform.unique.tsv"));
    assert_eq!(
        mappings,
        [
            "read_minus\tref_minus",
            "read_plus\tref_plus",
            "read_unknown\tref_unknown",
        ]
    );
    assert!(bed_names(&output_root.join("matrix_unused.bed")).is_empty());
    assert_eq!(
        bed_names(&output_root.join("matrix_novel.bed")),
        ["read_wrong_strand"]
    );

    let input_read_count = bed_names(&reads).len();
    let represented = mappings.len();
    let unused = bed_names(&output_root.join("matrix_unused.bed")).len();
    let novel = bed_names(&output_root.join("matrix_novel.bed")).len();
    assert_eq!(input_read_count, represented + unused + novel);
    assert!(
        (count_sum(&output_root.join("matrix_isoform_count.csv")) - represented as f64).abs()
            < 1e-9
    );
}

#[test]
fn flow_conserves_overlapping_gene_molecules_and_names_novel_isoforms_per_gene() {
    let fixture = TestDir::new("regression-overlapping-genes");
    let reads = fixture.path().join("reads.bed");
    let references = fixture.path().join("references.bed");
    let output_root = fixture.path().join("out");

    let reference_line = |name: &str, gene: &str| {
        format!(
            "chrF\t100\t250\t{name}\t100\t+\t100\t250\t0\t2\t50,50,\t0,100,\tnone\tnone\tnone\t-1,-1,\tisoform_anno\t{gene}\tgroup1\tnone\n"
        )
    };
    fs::write(
        &references,
        [
            reference_line("ref_a", "GENEA"),
            reference_line("ref_b", "GENEB"),
        ]
        .concat(),
    )
    .expect("write overlapping references");
    fs::write(
        &reads,
        [
            "chrF\t100\t250\tread_known\t1\t+\t100\t250\t0\t2\t50,50,\t0,100,\tnone\tnone\tnone\t-1,-1,\tnanopore_read\tnone\tgroup1\tnone\n",
            "chrF\t100\t250\tread_novel\t1\t+\t100\t250\t0\t2\t50,30,\t0,120,\tnone\tnone\tnone\t-1,-1,\tnanopore_read\tnone\tgroup1\tnone\n",
        ]
        .concat(),
    )
    .expect("write overlapping-gene reads");

    let output = run_trackcluster(&[
        "flow",
        "--reads",
        reads.to_str().unwrap(),
        "--reference",
        references.to_str().unwrap(),
        "--output-root",
        output_root.to_str().unwrap(),
        "--prefix",
        "overlap",
        "--threads",
        "1",
        "--junction-correction-min-support",
        "1",
        "--max-reads-per-gene",
        "0",
        "--force",
    ]);
    assert_success(&output, "overlapping-gene flow");

    for gene in ["GENEA", "GENEB"] {
        let gene_reads = trackcluster_rs::io::bed::read_bed12(
            output_root.join(gene).join(format!("{gene}_nano.bed")),
        )
        .unwrap()
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
        assert_eq!(gene_reads.len(), 2);
        assert!(gene_reads
            .iter()
            .all(|read| read.metadata().gene_id() == Some(gene)));
    }

    let catalog = trackcluster_rs::io::bed::read_bed12(output_root.join("overlap_isoform.bed"))
        .unwrap()
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    let novel = catalog
        .iter()
        .filter(|isoform| isoform.name.starts_with("tc_novel_v1:"))
        .collect::<Vec<_>>();
    assert_eq!(novel.len(), 2);
    assert_ne!(novel[0].name, novel[1].name);
    assert_eq!(
        novel
            .iter()
            .filter_map(|isoform| isoform.metadata().gene_id())
            .collect::<std::collections::BTreeSet<_>>(),
        std::collections::BTreeSet::from(["GENEA", "GENEB"])
    );

    let mappings = nonempty_lines(&output_root.join("overlap_read_to_isoform.unique.tsv"));
    let mut targets_by_read =
        std::collections::BTreeMap::<&str, std::collections::BTreeSet<&str>>::new();
    for line in &mappings {
        let (read, isoform) = line.split_once('\t').expect("mapping fields");
        targets_by_read.entry(read).or_default().insert(isoform);
    }
    assert_eq!(targets_by_read["read_known"], ["ref_a", "ref_b"].into());
    assert_eq!(targets_by_read["read_novel"].len(), 2);
    assert!(targets_by_read["read_novel"]
        .iter()
        .all(|isoform| isoform.starts_with("tc_novel_v1:")));

    let count_path = output_root.join("overlap_isoform_count.csv");
    let mut count_reader = csv::Reader::from_path(&count_path).unwrap();
    let counts = count_reader
        .records()
        .map(|record| {
            let record = record.unwrap();
            (record[1].to_owned(), record[2].parse::<f64>().unwrap())
        })
        .collect::<std::collections::BTreeMap<_, _>>();
    assert_eq!(counts.len(), 4);
    assert!(counts.values().all(|count| (*count - 0.5).abs() < 1e-9));
    assert!((count_sum(&count_path) - 2.0).abs() < 1e-9);
}

#[test]
fn desc_cli_reports_fusion_in_deterministic_gene_order() {
    let fixture = TestDir::new("regression-fusion");
    let isoforms = fixture.path().join("isoforms.bed");
    let references = fixture.path().join("references.bed");
    let prefix = fixture.path().join("annotation");

    fs::write(
        &references,
        [
            one_exon_bed_line("chr1", 210, 310, "ref_b", '+', "isoform_anno", "GENE_B"),
            one_exon_bed_line("chr1", 100, 200, "ref_a", '+', "isoform_anno", "GENE_A"),
        ]
        .concat(),
    )
    .expect("write fusion references");
    fs::write(
        &isoforms,
        one_exon_bed_line(
            "chr1",
            150,
            260,
            "fusion_read",
            '+',
            "nanopore_read",
            "none",
        ),
    )
    .expect("write fusion isoform");

    let output = run_trackcluster(&[
        "desc",
        "--isoform",
        isoforms.to_str().unwrap(),
        "--reference",
        references.to_str().unwrap(),
        "--out",
        prefix.to_str().unwrap(),
    ]);
    assert_success(&output, "fusion description");
    assert_eq!(
        nonempty_lines(&fixture.path().join("annotation_fusion.txt")),
        [
            "#schema\ttrackcluster-description-v2\tfusion",
            "isoform_id\tgene_ids",
            "fusion_read\tGENE_A;GENE_B",
        ]
    );
    assert!(
        nonempty_lines(&fixture.path().join("annotation_class12.txt"))
            .iter()
            .any(|line| line == "fusion_read\tfusion_gene")
    );
}

#[test]
fn desc_cli_uses_minus_reference_orientation_for_unknown_strand_isoforms() {
    let fixture = TestDir::new("regression-description-strands");
    let isoforms = fixture.path().join("isoforms.bed");
    let references = fixture.path().join("references.bed");
    let prefix = fixture.path().join("strand_annotation");

    fs::write(
        &references,
        reference_bed_line("chr1", 100, "ref_minus", '-', "GENE_MINUS"),
    )
    .expect("write minus-strand reference");
    fs::write(
        &isoforms,
        [
            one_exon_bed_line(
                "chr1",
                120,
                130,
                "query_unknown",
                '.',
                "nanopore_read",
                "GENE_MINUS",
            ),
            one_exon_bed_line(
                "chr1",
                120,
                130,
                "query_minus",
                '-',
                "nanopore_read",
                "GENE_MINUS",
            ),
        ]
        .concat(),
    )
    .expect("write minus and unknown queries");

    let output = run_trackcluster(&[
        "desc",
        "--isoform",
        isoforms.to_str().unwrap(),
        "--reference",
        references.to_str().unwrap(),
        "--out",
        prefix.to_str().unwrap(),
        "--offset-bp",
        "0",
    ]);
    assert_success(&output, "minus/unknown strand description");

    let rows = nonempty_lines(&fixture.path().join("strand_annotation_desc.txt"));
    assert_eq!(rows.len(), 4, "{rows:?}");
    let minus_fields = rows[2].split('\t').collect::<Vec<_>>();
    let unknown_fields = rows[3].split('\t').collect::<Vec<_>>();
    assert_eq!(minus_fields[0], "query_minus");
    assert_eq!(unknown_fields[0], "query_unknown");
    assert_eq!(minus_fields[1..], unknown_fields[1..]);
    assert!(
        rows[2].contains("5 primer miss") && rows[2].contains("3 primer miss"),
        "test must exercise a directional description rule: {:?}",
        rows[2]
    );
}

#[test]
fn flow_rejects_empty_reads_without_publishing_downstream_artifacts() {
    let fixture = TestDir::new("regression-empty");
    let reads = fixture.path().join("empty.bed");
    let references = fixture.path().join("references.bed");
    let output_root = fixture.path().join("out");
    fs::write(&reads, "").expect("write empty reads");
    fs::write(
        &references,
        reference_bed_line("chr1", 100, "ref", '+', "GENE"),
    )
    .expect("write reference");

    let output = run_trackcluster(&[
        "flow",
        "--reads",
        reads.to_str().unwrap(),
        "--reference",
        references.to_str().unwrap(),
        "--output-root",
        output_root.to_str().unwrap(),
        "--prefix",
        "empty",
        "--max-reads-per-gene",
        "0",
        "--force",
    ]);
    assert!(output.status.code().is_some_and(|code| code != 0));
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("no genes found"), "{stderr}");
    assert!(output_root.join("empty_gene.txt").exists());
    for downstream in [
        "empty_isoform.bed",
        "empty_unused.bed",
        "empty_read_to_isoform.tsv",
        "empty_read_to_isoform.unique.tsv",
        "empty_isoform_count.csv",
        "empty_desc.txt",
    ] {
        assert!(
            !output_root.join(downstream).exists(),
            "empty flow published downstream artifact {downstream}"
        );
    }
}
