#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 4 ]]; then
  echo "usage: $0 <reads.bed> <reference.bed> <output_root> <prefix> [threads] [sw_score] [batch_size] [batch_rounds]" >&2
  exit 2
fi

READS_BED="$1"
REFERENCE_BED="$2"
OUTPUT_ROOT="$3"
PREFIX="$4"
THREADS="${5:-30}"
SW_SCORE="${6:-11}"
BATCH_SIZE="${7:-2000}"
BATCH_ROUNDS="${8:-100}"

ts() { date +"%Y-%m-%d %H:%M:%S"; }

echo "[$(ts)] build: cargo build --release" >&2
cargo build --release

mkdir -p "$OUTPUT_ROOT"

LOG_FILE="$OUTPUT_ROOT/run.log"
# Mirror all output to a log file inside OUTPUT_ROOT.
exec > >(tee -a "$LOG_FILE") 2>&1

echo "[$(ts)] flow: output_root=$OUTPUT_ROOT prefix=$PREFIX threads=$THREADS sw_score=$SW_SCORE batch_size=$BATCH_SIZE batch_rounds=$BATCH_ROUNDS" >&2
target/release/trackcluster flow \
  --reads "$READS_BED" \
  --reference "$REFERENCE_BED" \
  --output-root "$OUTPUT_ROOT" \
  --prefix "$PREFIX" \
  --threads "$THREADS" \
  --sw-score "$SW_SCORE" \
  --batch-size "$BATCH_SIZE" \
  --batch-rounds "$BATCH_ROUNDS" \
  --force

echo "[$(ts)] done" >&2
