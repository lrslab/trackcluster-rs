use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::time::{SystemTime, UNIX_EPOCH};

use sha2::{Digest, Sha256};
use trackcluster_rs::model::Transcript;

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct StructureKey {
    chrom: String,
    strand: char,
    exons: Vec<(u32, u32)>,
}

impl StructureKey {
    fn from_transcript(transcript: &Transcript) -> Self {
        Self {
            chrom: transcript.chrom.clone(),
            strand: transcript.strand.as_char(),
            exons: transcript
                .exons
                .iter()
                .map(|exon| (exon.start.get(), exon.end.get()))
                .collect(),
        }
    }

    fn starts(&self) -> String {
        self.exons
            .iter()
            .map(|(start, _)| start.to_string())
            .collect::<Vec<_>>()
            .join(",")
    }

    fn ends(&self) -> String {
        self.exons
            .iter()
            .map(|(_, end)| end.to_string())
            .collect::<Vec<_>>()
            .join(",")
    }
}

struct TestDir(PathBuf);

impl TestDir {
    fn new(label: &str) -> Self {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("system clock before Unix epoch")
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "trackcluster-scientific-{label}-{}-{nonce}",
            std::process::id()
        ));
        fs::create_dir_all(&path).expect("create managed test directory");
        Self(path)
    }

    fn path(&self) -> &Path {
        &self.0
    }
}

impl Drop for TestDir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn repo_path(relative: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(relative)
}

fn read_transcripts(path: &Path) -> Vec<Transcript> {
    trackcluster_rs::io::bed::read_bed12(path)
        .unwrap_or_else(|error| panic!("open BED {path:?}: {error}"))
        .collect::<Result<Vec<_>, _>>()
        .unwrap_or_else(|error| panic!("parse BED {path:?}: {error}"))
}

fn run_success(command: &mut Command) -> Output {
    let output = command.output().expect("start trackcluster command");
    assert!(
        output.status.success(),
        "command failed with {}\nstdout:\n{}\nstderr:\n{}",
        output.status,
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn run_clusterj(reads: &Path, reference: &Path, output_bed: &Path, name2_mode: &str) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_trackcluster"));
    command
        .arg("clusterj")
        .arg("--reads")
        .arg(reads)
        .arg("--reference")
        .arg(reference)
        .arg("--out")
        .arg(output_bed)
        .arg("--junction-correction-min-support")
        .arg("1")
        .arg("--junction-correction-offset")
        .arg("0")
        .arg("--3prime-cluster-offset")
        .arg("0")
        .arg("--name2-mode")
        .arg(name2_mode);
    run_success(&mut command);
}

fn read_mapping(path: &Path) -> Vec<(String, String)> {
    let text = fs::read_to_string(path).expect("read mapping TSV");
    text.lines()
        .filter(|line| !line.trim().is_empty())
        .map(|line| {
            let (read_id, isoform_id) = line
                .split_once('\t')
                .unwrap_or_else(|| panic!("invalid mapping row {line:?}"));
            (read_id.to_owned(), isoform_id.to_owned())
        })
        .collect()
}

fn verify_sha256_manifest(directory: &Path) {
    let manifest_path = directory.join("SHA256SUMS");
    let manifest = fs::read_to_string(&manifest_path).expect("read SHA256SUMS");
    assert!(!manifest.trim().is_empty(), "empty {manifest_path:?}");

    for (line_index, line) in manifest.lines().enumerate() {
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let (expected, relative) = line.split_once("  ").unwrap_or_else(|| {
            panic!(
                "invalid SHA256SUMS row {} in {manifest_path:?}",
                line_index + 1
            )
        });
        assert_eq!(expected.len(), 64, "invalid SHA-256 width for {relative}");
        let bytes = fs::read(directory.join(relative))
            .unwrap_or_else(|error| panic!("read checksummed artifact {relative:?}: {error}"));
        let actual = format!("{:x}", Sha256::digest(bytes));
        assert_eq!(actual, expected, "checksum mismatch for {relative}");
    }
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct ParityRow {
    transcript_class: String,
    chrom: String,
    strand: char,
    exon_starts: String,
    exon_ends: String,
    represented_reads: String,
}

fn read_expected_parity(path: &Path) -> Vec<ParityRow> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .expect("open legacy expected projection");
    assert_eq!(
        reader
            .headers()
            .expect("read parity header")
            .iter()
            .collect::<Vec<_>>(),
        vec![
            "transcript_class",
            "chrom",
            "strand",
            "exon_starts_0based",
            "exon_ends_half_open",
            "represented_read_ids",
        ]
    );
    let mut rows = reader
        .records()
        .map(|record| {
            let record = record.expect("read expected parity record");
            ParityRow {
                transcript_class: record[0].to_owned(),
                chrom: record[1].to_owned(),
                strand: record[2].chars().next().expect("strand character"),
                exon_starts: record[3].to_owned(),
                exon_ends: record[4].to_owned(),
                represented_reads: record[5].to_owned(),
            }
        })
        .collect::<Vec<_>>();
    rows.sort();
    rows
}

#[test]
fn current_clusterj_matches_frozen_legacy_0_1_8_structures_and_membership() {
    let corpus = repo_path("tests/independent/legacy_trackcluster_0_1_8");
    verify_sha256_manifest(&corpus);

    let reads_path = corpus.join("inputs/reads.bed");
    let reference_path = corpus.join("inputs/reference.bed");
    let references = read_transcripts(&reference_path);
    let reference_structures = references
        .iter()
        .map(StructureKey::from_transcript)
        .collect::<BTreeSet<_>>();

    let output_dir = TestDir::new("legacy-parity");
    let output_bed = output_dir.path().join("isoforms.bed");
    run_clusterj(&reads_path, &reference_path, &output_bed, "full");

    let isoforms = read_transcripts(&output_bed);
    let mapping = read_mapping(&output_bed.with_extension("read_to_isoform.tsv"));
    let mut reads_by_isoform: HashMap<&str, BTreeSet<&str>> = HashMap::new();
    for (read_id, isoform_id) in &mapping {
        reads_by_isoform
            .entry(isoform_id)
            .or_default()
            .insert(read_id);
    }

    let mut actual = Vec::new();
    let mut structures_seen = BTreeSet::new();
    for isoform in &isoforms {
        let structure = StructureKey::from_transcript(isoform);
        assert!(
            structures_seen.insert(structure.clone()),
            "duplicate structural isoform in parity output: {structure:?}"
        );
        actual.push(ParityRow {
            transcript_class: if reference_structures.contains(&structure) {
                "known".to_owned()
            } else {
                "novel".to_owned()
            },
            chrom: structure.chrom.clone(),
            strand: structure.strand,
            exon_starts: structure.starts(),
            exon_ends: structure.ends(),
            represented_reads: reads_by_isoform
                .get(isoform.name.as_str())
                .map(|names| names.iter().copied().collect::<Vec<_>>().join(","))
                .unwrap_or_default(),
        });
    }
    actual.sort();

    let expected = read_expected_parity(&corpus.join("expected_structures.tsv"));
    assert_eq!(actual, expected);

    let input_read_ids = read_transcripts(&reads_path)
        .into_iter()
        .map(|read| read.name)
        .collect::<BTreeSet<_>>();
    let mapped_read_ids = mapping
        .iter()
        .map(|(read_id, _)| read_id.clone())
        .collect::<BTreeSet<_>>();
    assert_eq!(mapped_read_ids, input_read_ids);
    assert!(
        fs::read_to_string(output_bed.with_extension("unused.bed"))
            .expect("read unused BED")
            .is_empty(),
        "legacy parity corpus unexpectedly discarded a read"
    );
}

#[derive(Clone, Debug)]
struct TruthEntry {
    transcript_id: String,
    gene_id: String,
    transcript_class: String,
    expected_count: f64,
    support_bin: String,
    expression_bin: String,
    structure: StructureKey,
}

fn read_truth_entries(directory: &Path) -> Vec<TruthEntry> {
    let catalog = read_transcripts(&directory.join("truth_catalog.bed"));
    let structures_by_id = catalog
        .into_iter()
        .map(|transcript| {
            (
                transcript.name.clone(),
                StructureKey::from_transcript(&transcript),
            )
        })
        .collect::<HashMap<_, _>>();

    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(directory.join("truth_index.tsv"))
        .expect("open truth index");
    assert_eq!(
        reader
            .headers()
            .expect("read truth-index header")
            .iter()
            .collect::<Vec<_>>(),
        vec![
            "transcript_id",
            "gene_id",
            "transcript_class",
            "simulated_read_count",
            "support_bin",
            "expression_bin",
        ]
    );

    reader
        .records()
        .map(|record| {
            let record = record.expect("read truth-index record");
            let transcript_id = record[0].to_owned();
            let entry = TruthEntry {
                structure: structures_by_id
                    .get(&transcript_id)
                    .unwrap_or_else(|| panic!("truth ID {transcript_id:?} missing from catalog"))
                    .clone(),
                transcript_id,
                gene_id: record[1].to_owned(),
                transcript_class: record[2].to_owned(),
                expected_count: record[3].parse().expect("numeric simulated read count"),
                support_bin: record[4].to_owned(),
                expression_bin: record[5].to_owned(),
            };
            assert!(matches!(entry.transcript_class.as_str(), "known" | "novel"));
            entry
        })
        .collect()
}

fn read_counts(path: &Path) -> BTreeMap<String, (String, f64)> {
    let mut reader = csv::Reader::from_path(path).expect("open count CSV");
    assert_eq!(
        reader
            .headers()
            .expect("read count header")
            .iter()
            .collect::<Vec<_>>(),
        vec!["gene", "isoform_id", "count"]
    );
    let mut counts = BTreeMap::new();
    for record in reader.records() {
        let record = record.expect("read count record");
        let previous = counts.insert(
            record[1].to_owned(),
            (
                record[0].to_owned(),
                record[2].parse().expect("numeric isoform count"),
            ),
        );
        assert!(previous.is_none(), "duplicate count row for {}", &record[1]);
    }
    counts
}

#[derive(Clone, Debug)]
struct Prediction {
    isoform_id: String,
    gene_id: String,
    transcript_class: String,
    support_bin: String,
    expression_bin: String,
    observed_count: f64,
    structure: StructureKey,
}

#[derive(Clone, Copy)]
struct GroupSpec {
    scope: &'static str,
    transcript_class: &'static str,
    support_bin: &'static str,
    expression_bin: &'static str,
}

impl GroupSpec {
    fn includes(&self, transcript_class: &str, support_bin: &str, expression_bin: &str) -> bool {
        (self.transcript_class == "all" || self.transcript_class == transcript_class)
            && (self.support_bin == "all" || self.support_bin == support_bin)
            && (self.expression_bin == "all" || self.expression_bin == expression_bin)
    }
}

const GROUPS: &[GroupSpec] = &[
    GroupSpec {
        scope: "overall",
        transcript_class: "all",
        support_bin: "all",
        expression_bin: "all",
    },
    GroupSpec {
        scope: "class",
        transcript_class: "known",
        support_bin: "all",
        expression_bin: "all",
    },
    GroupSpec {
        scope: "class",
        transcript_class: "novel",
        support_bin: "all",
        expression_bin: "all",
    },
    GroupSpec {
        scope: "support",
        transcript_class: "all",
        support_bin: "high",
        expression_bin: "all",
    },
    GroupSpec {
        scope: "support",
        transcript_class: "all",
        support_bin: "low",
        expression_bin: "all",
    },
    GroupSpec {
        scope: "expression",
        transcript_class: "all",
        support_bin: "all",
        expression_bin: "high",
    },
    GroupSpec {
        scope: "expression",
        transcript_class: "all",
        support_bin: "all",
        expression_bin: "low",
    },
    GroupSpec {
        scope: "class_support_expression",
        transcript_class: "known",
        support_bin: "high",
        expression_bin: "high",
    },
    GroupSpec {
        scope: "class_support_expression",
        transcript_class: "known",
        support_bin: "low",
        expression_bin: "low",
    },
    GroupSpec {
        scope: "class_support_expression",
        transcript_class: "novel",
        support_bin: "high",
        expression_bin: "high",
    },
    GroupSpec {
        scope: "class_support_expression",
        transcript_class: "novel",
        support_bin: "low",
        expression_bin: "low",
    },
];

fn ratio(numerator: usize, denominator: usize) -> f64 {
    if denominator == 0 {
        0.0
    } else {
        numerator as f64 / denominator as f64
    }
}

fn render_metrics(truth: &[TruthEntry], predictions: &[Prediction]) -> String {
    let mut report = String::from(
        "scope\ttranscript_class\tsupport_bin\texpression_bin\ttruth_transcripts\tpredicted_transcripts\ttp\tfp\tfn\tprecision\trecall\tf1\texpected_abundance\tobserved_abundance\tmean_absolute_error\n",
    );

    for group in GROUPS {
        let truth_group = truth
            .iter()
            .filter(|entry| {
                group.includes(
                    &entry.transcript_class,
                    &entry.support_bin,
                    &entry.expression_bin,
                )
            })
            .collect::<Vec<_>>();
        let prediction_group = predictions
            .iter()
            .filter(|entry| {
                group.includes(
                    &entry.transcript_class,
                    &entry.support_bin,
                    &entry.expression_bin,
                )
            })
            .collect::<Vec<_>>();

        let truth_keys = truth_group
            .iter()
            .map(|entry| entry.structure.clone())
            .collect::<BTreeSet<_>>();
        let prediction_keys = prediction_group
            .iter()
            .map(|entry| entry.structure.clone())
            .collect::<BTreeSet<_>>();
        assert_eq!(
            truth_keys.len(),
            truth_group.len(),
            "duplicate truth structure"
        );
        assert_eq!(
            prediction_keys.len(),
            prediction_group.len(),
            "duplicate predicted structure"
        );

        let true_positives = truth_keys.intersection(&prediction_keys).count();
        let false_positives = prediction_keys.difference(&truth_keys).count();
        let false_negatives = truth_keys.difference(&prediction_keys).count();
        let precision = ratio(true_positives, true_positives + false_positives);
        let recall = ratio(true_positives, true_positives + false_negatives);
        let f1 = if precision + recall == 0.0 {
            0.0
        } else {
            2.0 * precision * recall / (precision + recall)
        };

        let expected_by_structure = truth_group
            .iter()
            .map(|entry| (entry.structure.clone(), entry.expected_count))
            .collect::<BTreeMap<_, _>>();
        let observed_by_structure = prediction_group
            .iter()
            .map(|entry| (entry.structure.clone(), entry.observed_count))
            .collect::<BTreeMap<_, _>>();
        let union = expected_by_structure
            .keys()
            .chain(observed_by_structure.keys())
            .cloned()
            .collect::<BTreeSet<_>>();
        let absolute_error_sum = union
            .iter()
            .map(|key| {
                (expected_by_structure.get(key).copied().unwrap_or(0.0)
                    - observed_by_structure.get(key).copied().unwrap_or(0.0))
                .abs()
            })
            .sum::<f64>();
        let mean_absolute_error = if union.is_empty() {
            0.0
        } else {
            absolute_error_sum / union.len() as f64
        };
        let expected_abundance = expected_by_structure.values().sum::<f64>();
        let observed_abundance = observed_by_structure.values().sum::<f64>();

        report.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{precision:.6}\t{recall:.6}\t{f1:.6}\t{expected_abundance:.6}\t{observed_abundance:.6}\t{mean_absolute_error:.6}\n",
            group.scope,
            group.transcript_class,
            group.support_bin,
            group.expression_bin,
            truth_group.len(),
            prediction_group.len(),
            true_positives,
            false_positives,
            false_negatives,
        ));
    }
    report
}

fn render_truth_keyed_tpm(truth: &[TruthEntry], predictions: &[Prediction]) -> String {
    let observed = predictions
        .iter()
        .map(|entry| (entry.structure.clone(), entry.observed_count))
        .collect::<BTreeMap<_, _>>();
    let total = observed.values().sum::<f64>();
    assert!(total > 0.0, "cannot calculate TPM from zero observations");

    let mut output = String::from("ID\tsimulated_rep1\n");
    for entry in truth {
        let count = observed.get(&entry.structure).copied().unwrap_or(0.0);
        output.push_str(&format!(
            "{}\t{:.6}\n",
            entry.transcript_id,
            count / total * 1_000_000.0
        ));
    }
    output
}

fn write_and_validate_lrgasp_artifacts(
    directory: &Path,
    isoforms: &[Transcript],
    mapping: &[(String, String)],
    predictions: &[Prediction],
) {
    let model_ids = isoforms
        .iter()
        .map(|transcript| transcript.name.as_str())
        .collect::<BTreeSet<_>>();

    let gtf = fs::read_to_string(directory.join("models.gtf")).expect("read exported GTF");
    let mut exon_model_ids = BTreeSet::new();
    let mut exon_rows = 0usize;
    for line in gtf.lines().filter(|line| !line.starts_with('#')) {
        let fields = line.split('\t').collect::<Vec<_>>();
        assert_eq!(fields.len(), 9, "invalid GTF row {line:?}");
        if fields[2] != "exon" {
            continue;
        }
        exon_rows += 1;
        // Attribute values are quoted and generated stable IDs may themselves
        // contain printable punctuation such as semicolons. Locate the quoted
        // value rather than naively splitting the attribute string on `;`.
        let transcript_id = fields[8]
            .split_once("transcript_id \"")
            .and_then(|(_, value)| value.split_once('"').map(|(id, _)| id))
            .unwrap_or_else(|| panic!("GTF exon lacks transcript_id: {line:?}"));
        assert!(
            fields[8].contains("gene_id \""),
            "GTF exon lacks gene_id: {line:?}"
        );
        exon_model_ids.insert(transcript_id);
    }
    assert_eq!(
        exon_rows,
        isoforms.iter().map(|tx| tx.exons.len()).sum::<usize>()
    );
    assert_eq!(exon_model_ids, model_ids);

    let read_model_map_path = directory.join("read_model_map.tsv");
    let mut map_writer = BufWriter::new(
        fs::File::create(&read_model_map_path).expect("create LRGASP read-model map"),
    );
    writeln!(map_writer, "read_id\ttranscript_id").expect("write map header");
    for (read_id, transcript_id) in mapping {
        assert!(model_ids.contains(transcript_id.as_str()));
        writeln!(map_writer, "{read_id}\t{transcript_id}").expect("write map row");
    }
    map_writer.flush().expect("flush LRGASP read-model map");
    let map_text = fs::read_to_string(&read_model_map_path).expect("read LRGASP map");
    assert_eq!(map_text.lines().next(), Some("read_id\ttranscript_id"));
    assert_eq!(map_text.lines().count(), mapping.len() + 1);

    let total_count = predictions
        .iter()
        .map(|prediction| prediction.observed_count)
        .sum::<f64>();
    let expression_path = directory.join("expression.tsv");
    let mut expression_writer = BufWriter::new(
        fs::File::create(&expression_path).expect("create LRGASP expression matrix"),
    );
    writeln!(expression_writer, "ID\tsimulated_rep1").expect("write expression header");
    for prediction in predictions {
        writeln!(
            expression_writer,
            "{}\t{:.6}",
            prediction.isoform_id,
            prediction.observed_count / total_count * 1_000_000.0
        )
        .expect("write expression row");
    }
    expression_writer
        .flush()
        .expect("flush LRGASP expression matrix");
    let mut expression_reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(&expression_path)
        .expect("open LRGASP expression matrix");
    assert_eq!(
        expression_reader
            .headers()
            .expect("read expression header")
            .iter()
            .collect::<Vec<_>>(),
        vec!["ID", "simulated_rep1"]
    );
    let mut tpm_sum = 0.0;
    for record in expression_reader.records() {
        let record = record.expect("read expression row");
        assert!(model_ids.contains(&record[0]));
        tpm_sum += record[1].parse::<f64>().expect("numeric TPM");
    }
    assert!((tpm_sum - 1_000_000.0).abs() < 0.01, "TPM sum={tpm_sum}");

    let sqanti = fs::read_to_string(directory.join("sqanti_input.tsv"))
        .expect("read SQANTI3 input audit table");
    let mut lines = sqanti.lines();
    assert_eq!(lines.next(), Some("#schema\ttrackcluster-sqanti-input-v1"));
    assert_eq!(
        lines.next(),
        Some("isoform_id\tgene_id\tchrom\tstrand\tlength\texon_count")
    );
    let sqanti_ids = lines
        .map(|line| line.split('\t').next().expect("SQANTI3 isoform ID"))
        .collect::<BTreeSet<_>>();
    assert_eq!(sqanti_ids, model_ids);
}

#[test]
fn simulated_truth_reports_known_novel_stratified_accuracy_and_standard_exports() {
    let corpus = repo_path("tests/independent/simulated_truth_v1");
    verify_sha256_manifest(&corpus);
    let truth = read_truth_entries(&corpus);
    let truth_structures = truth
        .iter()
        .map(|entry| entry.structure.clone())
        .collect::<BTreeSet<_>>();
    assert_eq!(truth_structures.len(), truth.len());

    let reference_path = corpus.join("reference.bed");
    let reads_path = corpus.join("reads.bed");
    let reference_structures = read_transcripts(&reference_path)
        .iter()
        .map(StructureKey::from_transcript)
        .collect::<BTreeSet<_>>();

    let output_dir = TestDir::new("simulated-truth");
    let output_bed = output_dir.path().join("isoforms.bed");
    run_clusterj(&reads_path, &reference_path, &output_bed, "coverage");
    let mapping_path = output_bed.with_extension("read_to_isoform.tsv");
    let mapping = read_mapping(&mapping_path);
    assert_eq!(
        mapping
            .iter()
            .map(|(read_id, _)| read_id)
            .collect::<BTreeSet<_>>()
            .len(),
        read_transcripts(&reads_path).len(),
        "each simulated molecule must be represented"
    );
    assert!(fs::read_to_string(output_bed.with_extension("unused.bed"))
        .expect("read simulated unused BED")
        .is_empty());

    let count_path = output_dir.path().join("counts.csv");
    let mut count_command = Command::new(env!("CARGO_BIN_EXE_trackcluster"));
    count_command
        .arg("count")
        .arg("--reads")
        .arg(&reads_path)
        .arg("--reference")
        .arg(&reference_path)
        .arg("--isoform")
        .arg(&output_bed)
        .arg("--read-to-isoform")
        .arg(&mapping_path)
        .arg("--assignment-mode")
        .arg("unique")
        .arg("--unique-assignment-junction-offset")
        .arg("0")
        .arg("--out")
        .arg(&count_path);
    run_success(&mut count_command);

    let isoforms = read_transcripts(&output_bed);
    let counts = read_counts(&count_path);
    assert_eq!(counts.len(), isoforms.len());
    let mut predictions = Vec::new();
    let mut prediction_structures = BTreeSet::new();
    for isoform in &isoforms {
        let structure = StructureKey::from_transcript(isoform);
        assert!(
            prediction_structures.insert(structure.clone()),
            "duplicate predicted structure: {structure:?}"
        );
        let (gene_id, observed_count) = counts
            .get(&isoform.name)
            .unwrap_or_else(|| panic!("missing count for {}", isoform.name));
        let support_bin = if *observed_count >= 5.0 {
            "high"
        } else {
            "low"
        };
        predictions.push(Prediction {
            isoform_id: isoform.name.clone(),
            gene_id: gene_id.clone(),
            transcript_class: if reference_structures.contains(&structure) {
                "known".to_owned()
            } else {
                "novel".to_owned()
            },
            support_bin: support_bin.to_owned(),
            expression_bin: support_bin.to_owned(),
            observed_count: *observed_count,
            structure,
        });
    }
    assert_eq!(prediction_structures, truth_structures);
    for prediction in &predictions {
        let truth_entry = truth
            .iter()
            .find(|entry| entry.structure == prediction.structure)
            .expect("prediction absent from truth catalog");
        assert_eq!(prediction.gene_id, truth_entry.gene_id);
    }

    let actual_metrics = render_metrics(&truth, &predictions);
    let expected_metrics = fs::read_to_string(corpus.join("expected_metrics.tsv"))
        .expect("read expected accuracy report");
    assert_eq!(actual_metrics, expected_metrics);

    let actual_tpm = render_truth_keyed_tpm(&truth, &predictions);
    let expected_tpm = fs::read_to_string(corpus.join("expected_expression_tpm.tsv"))
        .expect("read expected TPM report");
    assert_eq!(actual_tpm, expected_tpm);

    let models_gtf = output_dir.path().join("models.gtf");
    let sqanti_input = output_dir.path().join("sqanti_input.tsv");
    let mut export_command = Command::new(env!("CARGO_BIN_EXE_trackcluster"));
    export_command
        .arg("export")
        .arg("--input")
        .arg(&output_bed)
        .arg("--gtf")
        .arg(&models_gtf)
        .arg("--sqanti-input")
        .arg(&sqanti_input);
    run_success(&mut export_command);
    write_and_validate_lrgasp_artifacts(output_dir.path(), &isoforms, &mapping, &predictions);
}
