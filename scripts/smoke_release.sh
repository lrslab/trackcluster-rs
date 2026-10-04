#!/usr/bin/env bash
# Exercise the shipped examples with already-built binaries; no Cargo required.
set -euo pipefail

if [[ $# -ne 4 ]]; then
  echo "usage: $0 <bin-dir> <asset-root> <version> <work-dir>" >&2
  exit 2
fi

bin_dir="$(cd "$1" && pwd)"
asset_root="$(cd "$2" && pwd)"
version="$3"
mkdir -p "$4"
work_dir="$(cd "$4" && pwd)"
trackcluster="$bin_dir/trackcluster"

test "$("$trackcluster" --version)" = "trackcluster $version"
test "$("$bin_dir/clusterj_batch" --version)" = "clusterj-batch $version"
"$trackcluster" --help > "$work_dir/trackcluster-help.txt"
"$bin_dir/clusterj_batch" --help > "$work_dir/clusterj-batch-help.txt"
grep -F 'polya-aggregate' "$work_dir/trackcluster-help.txt"

"$trackcluster" validate-bed --input "$asset_root/examples/minimal.bed"
"$trackcluster" gff2bigg \
  --gff "$asset_root/examples/annotation.gff3" --out "$work_dir/reference.bed"
"$trackcluster" validate-bed --input "$work_dir/reference.bed"

"$trackcluster" flow \
  --reads "$asset_root/examples/reads.bed" --reference "$asset_root/examples/ref.bed" \
  --output-root "$work_dir/single" --prefix sample --threads 1 --heartbeat-seconds 0
"$trackcluster" flow \
  --manifest "$asset_root/examples/samples.tsv" --reference "$asset_root/examples/ref.bed" \
  --output-root "$work_dir/pooled" --prefix pooled --threads 1 --heartbeat-seconds 0
"$trackcluster" export --input "$work_dir/single/sample_isoform.bed" \
  --gtf "$work_dir/single/sample_isoform.gtf" --gff3 "$work_dir/single/sample_isoform.gff3"
"$trackcluster" flow --count-only \
  --reference "$asset_root/examples/ref.bed" --output-root "$work_dir/single" \
  --prefix sample --threads 1 --heartbeat-seconds 0

"$trackcluster" count \
  --reads "$asset_root/examples/reads.bed" --isoform "$asset_root/examples/ref.bed" \
  --assign-against-catalog --out "$work_dir/quant.csv"
"$trackcluster" count-multi \
  --manifest "$asset_root/examples/samples.tsv" --isoform "$asset_root/examples/ref.bed" \
  --assign-against-catalog --out "$work_dir/quant_multi"
"$trackcluster" polya-aggregate \
  --nanopolish "$asset_root/examples/polya.nanopolish.tsv" --sample S1 \
  --isoforms "$asset_root/examples/ref.bed" \
  --read-to-isoform "$work_dir/quant.read_to_isoform.tsv" --out "$work_dir/quant"
"$trackcluster" flow \
  --reads "$asset_root/examples/reads.bed" --reference "$asset_root/examples/ref.bed" \
  --output-root "$work_dir/polya" --prefix sample --threads 1 --heartbeat-seconds 0 \
  --polya-nanopolish "$asset_root/examples/polya.nanopolish.tsv" --polya-sample S1

python3 - "$work_dir" <<'PY'
import csv
import sys
from pathlib import Path

root = Path(sys.argv[1])


def rows(name, delimiter="\t"):
    with (root / name).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter=delimiter))


assert rows("quant.csv", ",") == [
    {"gene": "GENEA", "isoform_id": "ref_a", "count": "1"},
    {"gene": "GENEA", "isoform_id": "ref_b", "count": "0"},
]
assert rows("quant_multi.isoform_counts.matrix.tsv") == [
    {"gene": "GENEA", "isoform_id": "ref_a", "S1": "1", "S2": "1"},
    {"gene": "GENEA", "isoform_id": "ref_b", "S1": "0", "S2": "0"},
]
assert (root / "quant.read_to_isoform.tsv").read_text() == "read_trunc\tref_a\n"
assert (root / "quant_multi.read_to_isoform.tsv").read_text() == (
    "S1::read_s1\tref_a\nS2::read_s2\tref_a\n"
)
summary = {row["isoform_id"]: row for row in rows("quant.isoform_polya.tsv")}
assert summary["ref_a"]["caller"] == "nanopolish"
assert summary["ref_a"]["polya_reads"] == "1"
assert summary["ref_a"]["polya_median_nt"] == "80.250000"
assert summary["ref_b"]["assigned_reads"] == "0"
assert summary["ref_b"]["polya_median_nt"] == "NA"
calls = rows("polya/sample.read_polya.tsv")
assert len(calls) == 1
assert calls[0]["read_id"] == "read_trunc"
assert calls[0]["polya_length_nt"] == "80.25"
assert calls[0]["status"] == "estimated"
mapping = (root / "polya/sample_read_to_isoform.unique.tsv").read_text()
assert mapping == f"read_trunc\t{calls[0]['isoform_id']}\n"
assert rows("polya/sample.polya_qc.tsv")[0]["read_join_rate"] == "1.000000"
for name in [
    "single/sample_isoform.gtf",
    "single/sample_isoform.gff3",
    "single/sample_isoform_count.csv",
    "pooled/pooled.isoform_counts.matrix.tsv",
]:
    assert (root / name).stat().st_size > 0, name
print("Release examples passed: flow, export, recount, fixed catalog and Nanopolish poly(A).")
PY
