# 488 Demo Walkthrough

This document is a concrete, end-to-end walkthrough for the 488 multi-sample dataset under:

- `/Users/li/data/shoudong_488`

It shows both:

- a smoke run (small grouped inputs) for quick validation
- a full run (all grouped inputs) for production outputs

## 1) Data layout

Base files:

- `/Users/li/data/shoudong_488/reference.bed`
- `/Users/li/data/shoudong_488/reads.bed` (large pooled reads file)

Smoke inputs:

- `/Users/li/data/shoudong_488/multigroup_smoke_20260213_202018/samples.tsv`
- `/Users/li/data/shoudong_488/multigroup_smoke_20260213_202018/groups/*.bed`

Full inputs:

- `/Users/li/data/shoudong_488/multigroup_full_20260213_202410/samples.tsv`
- `/Users/li/data/shoudong_488/multigroup_full_20260213_202410/groups/*.bed`

## 2) Build binaries

From the repo root:

```bash
cargo build --release
```

You should have:

- `./target/release/trackcluster`
- `./target/release/clusterj_batch`

## 3) Optional: validate key BED inputs

```bash
./target/release/trackcluster validate-bed -i /Users/li/data/shoudong_488/reference.bed
./target/release/trackcluster validate-bed -i /Users/li/data/shoudong_488/multigroup_full_20260213_202410/groups/488_aba_1_s.bed
```

## 4) Smoke run (fast sanity check)

```bash
./target/release/trackcluster flow \
  --manifest /Users/li/data/shoudong_488/multigroup_smoke_20260213_202018/samples.tsv \
  --reference /Users/li/data/shoudong_488/reference.bed \
  --output-root /Users/li/data/shoudong_488/multigroup_smoke_20260213_202018/out \
  --prefix shoudong488 \
  --threads 8 \
  --sw-score -1
```

Expected key outputs:

- `/Users/li/data/shoudong_488/multigroup_smoke_20260213_202018/out/shoudong488_isoform.bed`
- `/Users/li/data/shoudong_488/multigroup_smoke_20260213_202018/out/shoudong488_isoform_count.csv`
- `/Users/li/data/shoudong_488/multigroup_smoke_20260213_202018/out/shoudong488.isoform_usage.long.tsv`
- `/Users/li/data/shoudong_488/multigroup_smoke_20260213_202018/out/shoudong488.isoform_counts.matrix.tsv`
- `/Users/li/data/shoudong_488/multigroup_smoke_20260213_202018/out/shoudong488.isoform_usage.group.tsv`
- `/Users/li/data/shoudong_488/multigroup_smoke_20260213_202018/out/clusterj_batch_summary.txt`

## 5) Full run (all grouped inputs)

```bash
./target/release/trackcluster flow \
  --manifest /Users/li/data/shoudong_488/multigroup_full_20260213_202410/samples.tsv \
  --reference /Users/li/data/shoudong_488/reference.bed \
  --output-root /Users/li/data/shoudong_488/multigroup_full_20260213_202410/out_grouped_rerun_20260213_203004 \
  --prefix shoudong488_full \
  --threads 8 \
  --sw-score -1
```

For reruns into the same output directory, add `--force`:

```bash
./target/release/trackcluster flow \
  --manifest /Users/li/data/shoudong_488/multigroup_full_20260213_202410/samples.tsv \
  --reference /Users/li/data/shoudong_488/reference.bed \
  --output-root /Users/li/data/shoudong_488/multigroup_full_20260213_202410/out_grouped_rerun_20260213_203004 \
  --prefix shoudong488_full \
  --threads 8 \
  --sw-score -1 \
  --force
```

## 6) Inspect run summary and outputs

Check batch summary:

```bash
sed -n '1,40p' /Users/li/data/shoudong_488/multigroup_full_20260213_202410/out_grouped_rerun_20260213_203004/clusterj_batch_summary.txt
```

Current recorded full-run summary includes:

- `total_genes: 25415`
- `processed: 25415`
- `errors: 0`
- `elapsed_seconds: 531.924451458` (with `threads: 8`)

Inspect output headers:

```bash
sed -n '1,5p' /Users/li/data/shoudong_488/multigroup_full_20260213_202410/out_grouped_rerun_20260213_203004/shoudong488_full.isoform_usage.long.tsv
sed -n '1,5p' /Users/li/data/shoudong_488/multigroup_full_20260213_202410/out_grouped_rerun_20260213_203004/shoudong488_full.isoform_counts.matrix.tsv
sed -n '1,5p' /Users/li/data/shoudong_488/multigroup_full_20260213_202410/out_grouped_rerun_20260213_203004/shoudong488_full.isoform_usage.group.tsv
sed -n '1,5p' /Users/li/data/shoudong_488/multigroup_full_20260213_202410/out_grouped_rerun_20260213_203004/shoudong488_full_isoform_count.csv
```

## 7) Manifest format used by this demo

`samples.tsv` (full) is:

```tsv
sample	group	reads
488_aba_1_s	488_aba	/Users/li/data/shoudong_488/multigroup_full_20260213_202410/groups/488_aba_1_s.bed
488_aba_2_s	488_aba	/Users/li/data/shoudong_488/multigroup_full_20260213_202410/groups/488_aba_2_s.bed
488_control_1_s	488_control	/Users/li/data/shoudong_488/multigroup_full_20260213_202410/groups/488_control_1_s.bed
488_control_2_s	488_control	/Users/li/data/shoudong_488/multigroup_full_20260213_202410/groups/488_control_2_s.bed
c24_aba_1_s	c24_aba	/Users/li/data/shoudong_488/multigroup_full_20260213_202410/groups/c24_aba_1_s.bed
c24_aba_2_s	c24_aba	/Users/li/data/shoudong_488/multigroup_full_20260213_202410/groups/c24_aba_2_s.bed
c24_control_1_s	c24_control	/Users/li/data/shoudong_488/multigroup_full_20260213_202410/groups/c24_control_1_s.bed
c24_control_2_s	c24_control	/Users/li/data/shoudong_488/multigroup_full_20260213_202410/groups/c24_control_2_s.bed
```

## 8) Notes

- `flow --manifest` performs pooled isoform discovery, then writes per-sample and per-group usage tables.
- Output file names are `<prefix>_*` and `<prefix>.*`, so keep `--prefix` stable if you want deterministic downstream paths.
- On Linux, prefer the `x86_64-unknown-linux-musl` release artifact to avoid glibc version mismatch issues.
