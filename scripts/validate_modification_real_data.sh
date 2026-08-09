#!/usr/bin/env bash
#
# Reproducible, opt-in validation against pinned public modification data.
#
# Usage:
#   scripts/validate_modification_real_data.sh m6anet
#   scripts/validate_modification_real_data.sh ont
#   scripts/validate_modification_real_data.sh all
#
# Environment:
#   MOD_VALIDATION_CACHE       download cache (default: target/real-data-cache)
#   MOD_VALIDATION_OUTPUT      generated outputs (default: target/real-data-validation)
#   TRACKCLUSTER_BIN           prebuilt binary; otherwise target/debug/trackcluster is built
#   M6ANET_REFERENCE_MODE      synthetic (default), ensembl91, gencode39, or provided
#   M6ANET_REFERENCE_GTF       existing GTF required by provided mode
#   ONT_RECORD_LIMIT           records streamed per public BAM (default: 1000)

set -euo pipefail

requested_suite="${1:-all}"
case "$requested_suite" in
    m6anet | ont | all) ;;
    *)
        echo "usage: $0 [m6anet|ont|all]" >&2
        exit 2
        ;;
esac

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cache_root="${MOD_VALIDATION_CACHE:-$repo_root/target/real-data-cache}"
output_root="${MOD_VALIDATION_OUTPUT:-$repo_root/target/real-data-validation}"
mkdir -p "$cache_root" "$output_root"
cd "$repo_root"

for required_command in curl python3; do
    if ! command -v "$required_command" >/dev/null 2>&1; then
        echo "required command not found: $required_command" >&2
        exit 1
    fi
done

if [[ -n "${TRACKCLUSTER_BIN:-}" ]]; then
    trackcluster_bin="$TRACKCLUSTER_BIN"
else
    cargo build --quiet --bin trackcluster
    trackcluster_bin="$repo_root/target/debug/trackcluster"
fi

sha256_file() {
    python3 - "$1" <<'PY'
import hashlib
import pathlib
import sys

source = pathlib.Path(sys.argv[1])
digest = hashlib.sha256()
with source.open("rb") as handle:
    for block in iter(lambda: handle.read(1024 * 1024), b""):
        digest.update(block)
print(digest.hexdigest())
PY
}

download_sha256() {
    local source_url="$1"
    local destination="$2"
    local expected_sha256="$3"
    local actual_sha256

    mkdir -p "$(dirname "$destination")"
    if [[ -f "$destination" ]]; then
        actual_sha256="$(sha256_file "$destination")"
    else
        actual_sha256=""
    fi
    if [[ "$actual_sha256" != "$expected_sha256" ]]; then
        curl --fail --location --retry 3 --connect-timeout 20 \
            --output "${destination}.part" "$source_url"
        mv "${destination}.part" "$destination"
        actual_sha256="$(sha256_file "$destination")"
    fi
    if [[ "$actual_sha256" != "$expected_sha256" ]]; then
        echo "SHA-256 mismatch for $destination" >&2
        echo "expected $expected_sha256" >&2
        echo "actual   $actual_sha256" >&2
        exit 1
    fi
}

verify_ensembl_bsd_sum() {
    python3 - "$1" <<'PY'
import pathlib
import sys

source = pathlib.Path(sys.argv[1])
checksum = 0
size = 0
with source.open("rb") as handle:
    for block in iter(lambda: handle.read(1024 * 1024), b""):
        size += len(block)
        for value in block:
            checksum = ((checksum >> 1) | ((checksum & 1) << 15))
            checksum = (checksum + value) & 0xffff
blocks = (size + 1023) // 1024
expected = (5485, 40882)
actual = (checksum, blocks)
if actual != expected:
    raise SystemExit(f"Ensembl CHECKSUMS mismatch: expected {expected}, got {actual}")
PY
}

run_m6anet_validation() {
    local upstream_commit="590ec277cb48d61774f0872395099e466022e810"
    local suite_cache="$cache_root/m6anet-$upstream_commit"
    local suite_output="$output_root/m6anet-$upstream_commit"
    local upstream_base="https://raw.githubusercontent.com/GoekeLab/m6anet/$upstream_commit/m6anet/tests/data"
    local indiv="$suite_cache/data.indiv_proba.csv.gz"
    local data_info="$suite_cache/data.info"
    local site_proba="$suite_cache/data.site_proba.csv.gz"
    local read_map="$suite_output/read_index_to_read_id.tsv"
    local synthetic_gtf="$suite_output/synthetic_projection.gtf"
    local reference_mode="${M6ANET_REFERENCE_MODE:-synthetic}"
    local reference_gtf
    local output_prefix="$suite_output/import"

    mkdir -p "$suite_cache" "$suite_output"
    download_sha256 \
        "$upstream_base/data.indiv_proba.csv.gz" \
        "$indiv" \
        "50f086b5c78ab703e3089401227e9e856b3541e5d4e8cb58a089309da0e5d44f"
    download_sha256 \
        "$upstream_base/data.info" \
        "$data_info" \
        "ce39729443224003aeee7085e6fe900c0db8a3f173f6f69056114e62e9892b3e"
    download_sha256 \
        "$upstream_base/data.site_proba.csv.gz" \
        "$site_proba" \
        "7b4b0e63eb0770213246610a5d49ff0857433118c520e5e01a8588d6a5ae3eec"

    python3 - "$indiv" "$data_info" "$read_map" "$synthetic_gtf" <<'PY'
import csv
import gzip
import pathlib
import sys

indiv_path = pathlib.Path(sys.argv[1])
data_info_path = pathlib.Path(sys.argv[2])
read_map_path = pathlib.Path(sys.argv[3])
gtf_path = pathlib.Path(sys.argv[4])
read_indexes = set()
maximum_position = {}
with gzip.open(indiv_path, "rt", newline="") as handle:
    for row in csv.DictReader(handle):
        read_indexes.add(row["read_index"])
        transcript_id = row["transcript_id"]
        position = int(row["transcript_position"])
        maximum_position[transcript_id] = max(maximum_position.get(transcript_id, -1), position)

with data_info_path.open(newline="") as handle:
    for row in csv.DictReader(handle):
        transcript_id = row["transcript_id"]
        position = int(row["transcript_position"])
        maximum_position[transcript_id] = max(maximum_position.get(transcript_id, -1), position)

with read_map_path.open("w", newline="") as handle:
    handle.write("read_index\tread_id\n")
    for read_index in sorted(read_indexes):
        handle.write(f"{read_index}\tm6anet-read-{read_index}\n")

# This reference intentionally validates source schema, versioned transcript IDs,
# and projection arithmetic only. A real annotation must pass the exact-version
# audit below before its genomic coordinates are used.
with gtf_path.open("w", newline="") as handle:
    for ordinal, transcript_id in enumerate(sorted(maximum_position), start=1):
        chrom = f"fixture_chr_{ordinal:03d}"
        gene_id = f"fixture_gene_{ordinal:03d}"
        end1 = maximum_position[transcript_id] + 1
        attributes = f'gene_id "{gene_id}"; transcript_id "{transcript_id}";'
        handle.write(
            f"{chrom}\tm6anet_fixture\ttranscript\t1\t{end1}\t.\t+\t.\t{attributes}\n"
        )
        handle.write(
            f"{chrom}\tm6anet_fixture\texon\t1\t{end1}\t.\t+\t.\t{attributes}\n"
        )
PY

    case "$reference_mode" in
        synthetic)
            reference_gtf="$synthetic_gtf"
            ;;
        ensembl91)
            reference_gtf="$suite_cache/Homo_sapiens.GRCh38.91.gtf.gz"
            if [[ ! -f "$reference_gtf" ]]; then
                curl --fail --location --retry 3 --connect-timeout 20 \
                    --output "${reference_gtf}.part" \
                    "https://ftp.ensembl.org/pub/release-91/gtf/homo_sapiens/Homo_sapiens.GRCh38.91.gtf.gz"
                mv "${reference_gtf}.part" "$reference_gtf"
            fi
            verify_ensembl_bsd_sum "$reference_gtf"
            ;;
        gencode39)
            reference_gtf="$suite_cache/gencode.v39.annotation.gtf.gz"
            download_sha256 \
                "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gtf.gz" \
                "$reference_gtf" \
                "bcb44a66c1cf567c8ebc941a12d3fe9565710d8aee10b0833a6f5b47d63c5c3a"
            ;;
        provided)
            if [[ -z "${M6ANET_REFERENCE_GTF:-}" ]]; then
                echo "M6ANET_REFERENCE_GTF is required with M6ANET_REFERENCE_MODE=provided" >&2
                exit 2
            fi
            reference_gtf="$M6ANET_REFERENCE_GTF"
            ;;
        *)
            echo "M6ANET_REFERENCE_MODE must be synthetic, ensembl91, gencode39, or provided" >&2
            exit 2
            ;;
    esac

    python3 - \
        "$indiv" \
        "$data_info" \
        "$reference_gtf" \
        "$suite_output/reference_compatibility.${reference_mode}.tsv" \
        "$suite_output/reference_compatibility_details.${reference_mode}.tsv" \
        "$reference_mode" <<'PY'
import csv
import gzip
import pathlib
import re
import sys

indiv_path = pathlib.Path(sys.argv[1])
data_info_path = pathlib.Path(sys.argv[2])
reference_path = pathlib.Path(sys.argv[3])
summary_path = pathlib.Path(sys.argv[4])
details_path = pathlib.Path(sys.argv[5])
reference_mode = sys.argv[6]

source_ids = set()
with gzip.open(indiv_path, "rt", newline="") as handle:
    source_ids.update(row["transcript_id"] for row in csv.DictReader(handle))
with data_info_path.open(newline="") as handle:
    source_ids.update(row["transcript_id"] for row in csv.DictReader(handle))

with reference_path.open("rb") as handle:
    is_gzip = handle.read(2) == b"\x1f\x8b"
open_text = gzip.open if is_gzip else open
transcript_pattern = re.compile(r'(?:^|;\s*)transcript_id\s+"([^"]+)"')
annotation_ids = set()
with open_text(reference_path, "rt") as handle:
    for line_number, line in enumerate(handle, start=1):
        if not line or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9:
            raise SystemExit(
                f"reference audit: {reference_path}:{line_number} has {len(fields)} fields, expected 9"
            )
        if fields[2].lower() not in {"transcript", "exon"}:
            continue
        match = transcript_pattern.search(fields[8])
        if match:
            annotation_ids.add(match.group(1))

def base_id(transcript_id):
    base, separator, version = transcript_id.rpartition(".")
    return base if separator and base and version.isdigit() else transcript_id

versions_by_base = {}
for transcript_id in annotation_ids:
    versions_by_base.setdefault(base_id(transcript_id), []).append(transcript_id)
for versions in versions_by_base.values():
    versions.sort()

details = []
for source_id in sorted(source_ids):
    if source_id in annotation_ids:
        status = "exact"
        alternatives = source_id
    else:
        available = versions_by_base.get(base_id(source_id), [])
        status = "version_mismatch" if available else "absent"
        alternatives = ",".join(available) if available else "NA"
    details.append((source_id, status, alternatives))

counts = {
    "exact": sum(status == "exact" for _, status, _ in details),
    "version_mismatch": sum(status == "version_mismatch" for _, status, _ in details),
    "absent": sum(status == "absent" for _, status, _ in details),
}
with details_path.open("w", newline="") as handle:
    handle.write("source_transcript_id\tstatus\tannotation_transcript_ids\n")
    for row in details:
        handle.write("\t".join(map(str, row)) + "\n")
with summary_path.open("w", newline="") as handle:
    handle.write("metric\tvalue\n")
    handle.write(f"reference_mode\t{reference_mode}\n")
    handle.write(f"source_transcripts\t{len(source_ids)}\n")
    for status in ("exact", "version_mismatch", "absent"):
        handle.write(f"{status}\t{counts[status]}\n")

if counts["exact"] != len(source_ids):
    raise SystemExit(
        "m6Anet reference audit failed before projection: "
        f"{counts['exact']}/{len(source_ids)} exact, "
        f"{counts['version_mismatch']} version mismatches, {counts['absent']} absent; "
        f"see {summary_path} and {details_path}"
    )
PY

    "$trackcluster_bin" mod-import-m6anet \
        --sample official_fixture \
        --assay-id m6anet_hct116_rna002 \
        --indiv "$indiv" \
        --data-info "$data_info" \
        --site-proba "$site_proba" \
        --read-map "$read_map" \
        --reference "$reference_gtf" \
        --model-id HCT116_RNA002 \
        --caller-version 2.1.0 \
        --candidate-rule DRACH \
        --min-reads 20 \
        --input-format gtf \
        --out "$output_prefix"

    python3 - \
        "$output_prefix.import_qc.tsv" \
        "$output_prefix.observations.tsv" \
        "$suite_output/validation_report.tsv" \
        "$reference_mode" \
        "$upstream_commit" \
        "$(sha256_file "$reference_gtf")" <<'PY'
import csv
import pathlib
import sys

qc_path = pathlib.Path(sys.argv[1])
observations_path = pathlib.Path(sys.argv[2])
report_path = pathlib.Path(sys.argv[3])
reference_mode = sys.argv[4]
upstream_commit = sys.argv[5]
reference_sha256 = sys.argv[6]

with qc_path.open(newline="") as handle:
    qc = {row["metric"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")}
expected = {
    "input_rows": "5595",
    "unique_observations": "5595",
    "duplicate_exact": "0",
    "read_map_entries": "4245",
    "read_map_entries_used": "4245",
    "source_transcripts": "52",
    "projection_transcripts_loaded": "82",
    "source_sites": "101",
    "data_info_sites": "248",
    "data_info_retained_sites": "101",
    "data_info_filtered_sites": "147",
    "data_info_total_reads": "7000",
    "data_info_retained_reads": "5595",
    "data_info_minimum_reads": "20",
    "site_probability_sites": "101",
    "site_probability_total_reads": "5595",
    "site_probability_sites_at_or_above_threshold": "35",
}
for metric, expected_value in expected.items():
    actual_value = qc.get(metric)
    if actual_value != expected_value:
        raise SystemExit(
            f"m6Anet validation failed: {metric} expected {expected_value}, got {actual_value}"
        )

with observations_path.open(newline="") as handle:
    observation_rows = sum(1 for _ in csv.DictReader(handle, delimiter="\t"))
if observation_rows != 5595:
    raise SystemExit(
        f"m6Anet validation failed: expected 5595 observation rows, got {observation_rows}"
    )

report = [
    ("metric", "value"),
    ("status", "pass"),
    ("source", "GoekeLab/m6anet official test fixture"),
    ("upstream_commit", upstream_commit),
    ("projection_reference", reference_mode),
    ("projection_reference_sha256", reference_sha256),
]
report.extend((metric, value) for metric, value in expected.items())
with report_path.open("w", newline="") as handle:
    for metric, value in report:
        handle.write(f"{metric}\t{value}\n")
print(f"m6Anet real-data validation passed: {report_path}")
PY
}

extract_remote_bam_prefix() {
    local source_url="$1"
    local destination="$2"
    local record_limit="$3"
    local existing_count=""

    if [[ -f "$destination" ]] && samtools quickcheck "$destination" 2>/dev/null; then
        existing_count="$(samtools view -c "$destination")"
    fi
    if [[ "$existing_count" == "$record_limit" ]]; then
        return
    fi

    set +e
    set +o pipefail
    samtools view -h "$source_url" |
        awk -v limit="$record_limit" \
            'BEGIN { n = 0 } /^@/ { print; next } n < limit { print; n++ } n == limit { exit }' |
        samtools view -b -o "${destination}.part" -
    local pipeline_codes=("${PIPESTATUS[@]}")
    set -o pipefail
    set -e
    if [[ "${pipeline_codes[2]}" -ne 0 || "${pipeline_codes[1]}" -ne 0 ]] ||
        [[ "${pipeline_codes[0]}" -ne 0 && "${pipeline_codes[0]}" -ne 141 ]]; then
        echo "failed to stream BAM prefix from $source_url: ${pipeline_codes[*]}" >&2
        exit 1
    fi
    mv "${destination}.part" "$destination"
    samtools quickcheck "$destination"
    existing_count="$(samtools view -c "$destination")"
    if [[ "$existing_count" != "$record_limit" ]]; then
        echo "expected $record_limit BAM records, got $existing_count in $destination" >&2
        exit 1
    fi
}

run_ont_validation() {
    local release="rna-mod-validation-all5mer-2026.07"
    local source_root="https://ont-open-data.s3.amazonaws.com/$release"
    local suite_cache="$cache_root/$release"
    local suite_output="$output_root/$release"
    local record_limit="${ONT_RECORD_LIMIT:-1000}"

    if ! command -v samtools >/dev/null 2>&1; then
        echo "ONT validation requires samtools for bounded remote BAM streaming" >&2
        exit 1
    fi
    if [[ ! "$record_limit" =~ ^[1-9][0-9]*$ ]] || ((record_limit < 1000)); then
        echo "ONT_RECORD_LIMIT must be an integer of at least 1000" >&2
        exit 2
    fi

    mkdir -p "$suite_cache/references" "$suite_output"
    download_sha256 \
        "$source_root/README.txt" \
        "$suite_cache/README.txt" \
        "db427e2c5e7e158a0fc3dc2c7a12e07f9cfaeb4d1bb0a87a816e8616758cf872"
    download_sha256 \
        "$source_root/commands.sh" \
        "$suite_cache/commands.sh" \
        "1f86a9a98ad8deb74897200104c3f2128f4b00ba5700198dde43d34ab8fbde51"
    download_sha256 \
        "$source_root/references/A_reference.fasta" \
        "$suite_cache/references/A_reference.fasta" \
        "c1170a99fafadbd0aac6d64d7842588bbcdb8c163780b86aa71dfa836a161700"
    download_sha256 \
        "$source_root/references/control_A_positions.bed" \
        "$suite_cache/references/control_A_positions.bed" \
        "b1b2257d66f87ae1f0a89c833c99919c33a9b31a9c94b24ba3ebc1217f13beba"
    download_sha256 \
        "$source_root/references/m6A_positions.bed" \
        "$suite_cache/references/m6A_positions.bed" \
        "34b6387ca4b1f18fd0b86dcfcab143d9cdc436c3652bddd408516d82fae96e2a"

    for sample_name in control_A m6A; do
        local bam_prefix="$suite_output/${sample_name}.first${record_limit}.bam"
        local output_prefix="$suite_output/$sample_name"
        extract_remote_bam_prefix \
            "$source_root/basecalls/${sample_name}.bam" \
            "$bam_prefix" \
            "$record_limit"
        "$trackcluster_bin" mod-import-dorado \
            --sample "$sample_name" \
            --assay-id ont_all5mer_v6_m6a \
            --bam "$bam_prefix" \
            --mod-code A+a \
            --model-id rna004_sup_v6.0.0_inosine_m6A_2OmeA_v1 \
            --chemistry RNA004 \
            --caller-version 2.0.0+20e87c8b \
            --candidate-rule all-target-canonical-bases \
            --source-emission-threshold 0.05 \
            --out "$output_prefix"
    done

    python3 - \
        "$suite_output/control_A.observations.tsv" \
        "$suite_output/m6A.observations.tsv" \
        "$suite_output/control_A.import_qc.tsv" \
        "$suite_output/m6A.import_qc.tsv" \
        "$suite_cache/references/control_A_positions.bed" \
        "$suite_cache/references/m6A_positions.bed" \
        "$suite_output/validation_report.tsv" \
        "$record_limit" <<'PY'
import csv
import pathlib
import statistics
import sys

control_observations = pathlib.Path(sys.argv[1])
modified_observations = pathlib.Path(sys.argv[2])
control_qc_path = pathlib.Path(sys.argv[3])
modified_qc_path = pathlib.Path(sys.argv[4])
control_truth = pathlib.Path(sys.argv[5])
modified_truth = pathlib.Path(sys.argv[6])
report_path = pathlib.Path(sys.argv[7])
record_limit = int(sys.argv[8])

def read_qc(source):
    with source.open(newline="") as handle:
        return {row["metric"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")}

def truth_metrics(observations_path, truth_path, expected_modified):
    truth = set()
    with truth_path.open() as handle:
        for line in handle:
            chrom, start, end, *_ = line.rstrip("\n").split("\t")
            if int(end) != int(start) + 1:
                raise SystemExit(f"truth interval is not one base: {line.rstrip()}")
            truth.add((chrom, int(start)))
    by_site = {site: 0 for site in truth}
    hard_calls = []
    explicit_probabilities = []
    unknown = 0
    with observations_path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            site = (row["chrom"], int(row["pos0"]))
            if site not in truth:
                continue
            state = row["observation_state"]
            if state == "explicit_probability":
                probability = float(row["probability"])
                explicit_probabilities.append(probability)
                hard_call = probability >= 0.5
            elif state == "implicit_below_emission_threshold":
                hard_call = False
            elif state == "unknown":
                unknown += 1
                continue
            else:
                raise SystemExit(f"unexpected observation state {state!r}")
            hard_calls.append(hard_call)
            by_site[site] += 1
    if not hard_calls:
        raise SystemExit(f"no truth-site calls in {observations_path}")
    correct = sum(call == expected_modified for call in hard_calls)
    return {
        "truth_sites": len(truth),
        "covered_truth_sites": sum(count > 0 for count in by_site.values()),
        "truth_site_calls": len(hard_calls),
        "truth_site_unknown": unknown,
        "hard_call_accuracy": correct / len(hard_calls),
        "modified_fraction": sum(hard_calls) / len(hard_calls),
        "mean_explicit_probability": statistics.mean(explicit_probabilities),
        "median_site_coverage": statistics.median(by_site.values()),
    }

control_qc = read_qc(control_qc_path)
modified_qc = read_qc(modified_qc_path)
for label, qc in (("control_A", control_qc), ("m6A", modified_qc)):
    if qc.get("total_records") != str(record_limit):
        raise SystemExit(f"{label}: total_records does not match requested prefix")
    if qc.get("candidate_observations_complete") != "true":
        raise SystemExit(f"{label}: candidate universe was marked incomplete")
    if qc.get("unknown_candidates") != "0":
        raise SystemExit(f"{label}: unexpected unknown candidates")

control = truth_metrics(control_observations, control_truth, False)
modified = truth_metrics(modified_observations, modified_truth, True)
for label, metrics in (("control_A", control), ("m6A", modified)):
    if metrics["truth_sites"] != 256 or metrics["covered_truth_sites"] != 256:
        raise SystemExit(f"{label}: expected all 256 truth sites to be covered")
if control["hard_call_accuracy"] < 0.90:
    raise SystemExit(f"control_A accuracy below 0.90: {control['hard_call_accuracy']}")
if modified["hard_call_accuracy"] < 0.80:
    raise SystemExit(f"m6A accuracy below 0.80: {modified['hard_call_accuracy']}")
direction_delta = modified["modified_fraction"] - control["modified_fraction"]
if direction_delta < 0.70:
    raise SystemExit(f"m6A-control modified-fraction delta below 0.70: {direction_delta}")

rows = [
    ("metric", "value"),
    ("status", "pass"),
    ("source", "ONT RNA all-5-mer synthetic ground truth"),
    ("release", "rna-mod-validation-all5mer-2026.07"),
    ("records_per_bam", str(record_limit)),
    ("source_control_A_bam_bytes", "108095235"),
    ("source_control_A_bam_etag", "252abe4693cb05119e5b955fbceb0f80-13"),
    ("source_m6A_bam_bytes", "123340699"),
    ("source_m6A_bam_etag", "42e444584a3c3ae19df21efb8c1047fc-15"),
]
for label, metrics in (("control_A", control), ("m6A", modified)):
    for metric, value in metrics.items():
        rows.append((f"{label}_{metric}", str(value)))
rows.append(("m6A_minus_control_modified_fraction", str(direction_delta)))
with report_path.open("w", newline="") as handle:
    for metric, value in rows:
        handle.write(f"{metric}\t{value}\n")
print(f"ONT real-data validation passed: {report_path}")
PY
}

case "$requested_suite" in
    m6anet)
        run_m6anet_validation
        ;;
    ont)
        run_ont_validation
        ;;
    all)
        run_m6anet_validation
        run_ont_validation
        ;;
esac
