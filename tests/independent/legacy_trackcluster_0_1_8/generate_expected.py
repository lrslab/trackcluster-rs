#!/usr/bin/env python3
"""Regenerate the frozen TrackCluster 0.1.8 structural parity projection.

This script intentionally imports and executes the legacy implementation. It is
not used by CI; CI consumes only the checked-in projection. Set
TRACKCLUSTER_LEGACY_ROOT to a TrackCluster 0.1.8 source checkout.
"""

from __future__ import annotations

import csv
import hashlib
import os
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
LEGACY_ROOT = Path(os.environ["TRACKCLUSTER_LEGACY_ROOT"]).resolve()

for checksum_line in (HERE / "LEGACY_SOURCE_SHA256SUMS").read_text().splitlines():
    expected_checksum, relative_path = checksum_line.split("  ", 1)
    source_path = LEGACY_ROOT / relative_path
    actual_checksum = hashlib.sha256(source_path.read_bytes()).hexdigest()
    if actual_checksum != expected_checksum:
        raise SystemExit(
            f"legacy source checksum mismatch for {relative_path}: "
            f"expected {expected_checksum}, found {actual_checksum}"
        )

sys.path.insert(0, str(LEGACY_ROOT))

import trackcluster  # noqa: E402
from trackcluster.clusterj import flow_junction_cluster  # noqa: E402
from trackcluster.tracklist import read_bigg  # noqa: E402


if trackcluster.__version__ != "0.1.8":
    raise SystemExit(
        f"expected TrackCluster 0.1.8, found {trackcluster.__version__!r}"
    )


references = read_bigg(str(HERE / "inputs" / "reference.bed"))
reads = read_bigg(str(HERE / "inputs" / "reads.bed"))
reference_names = {reference.name for reference in references}
result = flow_junction_cluster(reads, references)

rows = []
for transcript in result:
    transcript.get_exon()
    represented_reads = set(transcript.subread)
    if transcript.name not in reference_names:
        # Legacy name2/subread contains reads merged *into* the representative,
        # but not the representative read itself.
        represented_reads.add(transcript.name)
    rows.append(
        (
            "known" if transcript.name in reference_names else "novel",
            transcript.chrom,
            transcript.strand,
            ",".join(str(start) for start, _ in transcript.exon),
            ",".join(str(end) for _, end in transcript.exon),
            ",".join(sorted(represented_reads)),
        )
    )

rows.sort()
with (HERE / "expected_structures.tsv").open("w", newline="") as stream:
    writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
    writer.writerow(
        [
            "transcript_class",
            "chrom",
            "strand",
            "exon_starts_0based",
            "exon_ends_half_open",
            "represented_read_ids",
        ]
    )
    writer.writerows(rows)
