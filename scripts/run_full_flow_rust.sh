#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if [[ $# -lt 4 ]]; then
  echo "usage: $0 <reads.bed> <reference.bed> <output_root> <prefix> [threads] [sw_score] [batch_size] [batch_rounds] [name2_mode] [max_reads_per_gene] [heartbeat_seconds] [heartbeat_top]" >&2
  exit 2
fi

READS_BED="$1"
REFERENCE_BED="$2"
OUTPUT_ROOT="$3"
PREFIX="$4"
THREADS="${5:-8}"
SW_SCORE="${6:--1}"
BATCH_SIZE="${7:-500}"
BATCH_ROUNDS="${8:-100}"
NAME2_MODE="${9:-coverage}"
MAX_READS_PER_GENE="${10:-5000}"
HEARTBEAT_SECONDS="${11:-60}"
HEARTBEAT_TOP="${12:-5}"

ts() { date +"%Y-%m-%d %H:%M:%S"; }

if [[ -x "$root/trackcluster" ]]; then
  TRACKCLUSTER="$root/trackcluster"
  echo "[$(ts)] binary: $TRACKCLUSTER" >&2
else
  echo "[$(ts)] build: cargo build --release" >&2
  cargo build --release
  TRACKCLUSTER="$root/target/release/trackcluster"
fi

mkdir -p "$OUTPUT_ROOT"

LOG_FILE="$OUTPUT_ROOT/run.log"
# Mirror all output to a log file inside OUTPUT_ROOT.
exec > >(tee -a "$LOG_FILE") 2>&1

echo "[$(ts)] flow: output_root=$OUTPUT_ROOT prefix=$PREFIX threads=$THREADS sw_score=$SW_SCORE batch_size=$BATCH_SIZE batch_rounds=$BATCH_ROUNDS name2_mode=$NAME2_MODE max_reads_per_gene=$MAX_READS_PER_GENE heartbeat_seconds=$HEARTBEAT_SECONDS heartbeat_top=$HEARTBEAT_TOP" >&2
"$TRACKCLUSTER" flow \
  --reads "$READS_BED" \
  --reference "$REFERENCE_BED" \
  --output-root "$OUTPUT_ROOT" \
  --prefix "$PREFIX" \
  --threads "$THREADS" \
  --sw-score "$SW_SCORE" \
  --batch-size "$BATCH_SIZE" \
  --batch-rounds "$BATCH_ROUNDS" \
  --name2-mode "$NAME2_MODE" \
  --max-reads-per-gene "$MAX_READS_PER_GENE" \
  --heartbeat-seconds "$HEARTBEAT_SECONDS" \
  --heartbeat-top "$HEARTBEAT_TOP" \
  --force

echo "[$(ts)] done" >&2
