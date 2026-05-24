#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"

tmp_dir="tmp/golden-gen"
rm -rf "$tmp_dir"
mkdir -p "$tmp_dir"

cargo run -- clusterj \
  -s tests/fixtures/reads.bed \
  -r tests/fixtures/ref.bed \
  -o "$tmp_dir/isoform.bed"

cargo run -- count \
  -s tests/fixtures/reads.bed \
  -r tests/fixtures/ref.bed \
  -i "$tmp_dir/isoform.bed" \
  --out "$tmp_dir/isoform_count.csv"

mkdir -p tests/golden/clusterj tests/golden/count
cp "$tmp_dir/isoform.bed" tests/golden/clusterj/isoform.bed
cp "$tmp_dir/isoform.read_to_isoform.tsv" tests/golden/clusterj/isoform.read_to_isoform.tsv
cp "$tmp_dir/isoform.unused.bed" tests/golden/clusterj/isoform.unused.bed
cp "$tmp_dir/isoform_count.csv" tests/golden/count/isoform_count.csv

echo "Updated goldens in tests/golden/"
