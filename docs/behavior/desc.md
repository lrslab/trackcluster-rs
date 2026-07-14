# `desc` and Figure 2 classification contract

This document defines the Rust implementation of the TrackCluster paper's
Figure 2 classification. The output labels remain compatible with the legacy
Python implementation in `trackcluster/post.py` and `trackcluster/flow.py`.

## Inputs

- Isoform-like BED, usually produced by `clusterj`, `cluster`, or `addgene`.
- Reference BED containing `isoform_anno` records with gene names.

Each input transcript needs a populated gene-name field for reference
comparison. The Rust implementation reads and writes this value through the
TrackCluster bigGenePred extra-field convention.

## Outputs

The Rust `desc` command writes four `trackcluster-description-v2` files. The
prefix is controlled by `-o/--out`.

- `<prefix>_desc.txt`: `isoform_id`, `reference_id`, `gene_id`,
  `missing_features`, `extra_features`
- `<prefix>_class4.txt`: `isoform_id`, `class`
- `<prefix>_fusion.txt`: `isoform_id`, `gene_ids` (only multi-gene entries)
- `<prefix>_class12.txt`: `isoform_id`, `class`

## Reference selection

An isoform is compared only with references on the same chromosome and a
compatible strand. Its nearest reference minimizes these values in order:

1. number of extra boundary groups;
2. number of extra boundaries;
3. number of missed boundary groups;
4. number of missed boundaries;
5. total transcript-end displacement;
6. reference ID, for a deterministic final tie-break.

Boundary vectors are placed in transcript order. Consequently, 5′ and 3′
always mean biological transcript direction: their genomic sides are reversed
on the minus strand.

## Figure 2 `class12`

`class12` consists of the 11 novel-event labels in Figure 2 G–M plus the
compatibility bucket `reference`. These are the only valid values.

| Figure panel | Figure meaning | Stable `class12` value |
| --- | --- | --- |
| G | 3′ terminal exon addition | `3'extra` |
| G | 3′ terminal exon loss | `3'missing` |
| H | 5′ terminal exon addition | `5'extra` |
| H | 5′ terminal exon loss | `5'missing` |
| I | Extra UTR with a matched splice chain | `full_matched>=` |
| I | Missing UTR with a matched splice chain | `full_matched<` |
| J | Internal exon inclusion | `inner_extra_exon` |
| J | Internal exon skipping | `inner_miss_exon` |
| K | Intron retention | `intron_retention` |
| L | Alternative splice site/new junction | `new_junction` |
| M | Fusion of two or more genes | `fusion_gene` |
| — | Input record is annotated (`isoform_anno`) | `reference` |

The legacy strings are intentionally retained for downstream compatibility.
For example, `full_matched>=` means Figure 2 “Extra UTR”; it does not mean that
an exactly equal, non-reference transcript automatically qualifies as novel.

Complete terminal-exon addition/loss is distinct from moving the outer
boundary of the same terminal exon. A one-for-one replacement of splice-site
boundaries is Figure 2L `new_junction`, not exon inclusion/skipping or intron
retention. Intron retention is detected geometrically when one query exon spans
a reference intron, including the first or last intron.

### Figure 2I UTR threshold

Extra/Missing UTR requires all of the following:

- the internal splice-boundary chain matches the selected reference within
  `--offset-bp`;
- the difference is confined to transcript ends;
- the absolute difference in summed exon length is at least 5% of the selected
  reference's summed exon length.

The threshold is evaluated without floating-point rounding:

```text
abs(query_summed_exon_length - reference_summed_exon_length) * 100
    >= reference_summed_exon_length * 5
```

Exactly 5% qualifies. A longer query is `full_matched>=`; a shorter query is
`full_matched<`. An exactly equal query or a difference below 5% is not a novel
UTR event. Upstream clustering is expected to merge such sub-threshold tracks;
when standalone `desc` receives one, it has no `class12` row unless another
Figure 2 rule applies. An annotated input record still becomes `reference`.

Figure 2 does not separately specify a case in which one end extends while the
other truncates. The implementation uses the net difference in summed exon
length; if that net difference is below 5%, it does not emit a UTR bucket.

### Overlap resolution: later rules overwrite earlier rules

Evidence can match more than one description. `class12` evaluates the following
list from low to high priority:

```text
new_junction
< 5'missing
< 3'missing
< 5'extra
< 3'extra
< intron_retention
< inner_miss_exon
< inner_extra_exon
< fusion_gene
< full_matched<
< full_matched>=
< reference
```

The evaluator deliberately continues after a match. Every later match
overwrites the previous value, and the last matching label is written. For
example, an isoform missing complete exons at both transcript ends matches both
terminal-missing rules and is written as `3'missing`; `reference` always has
the highest priority. Figure-specific evidence is made disjoint where needed:
an alternative splice-site replacement is removed from whole-exon/retention
evidence, and a fusion is not also treated as a UTR event.

An isoform that matches no Figure 2 evidence set is omitted from
`<prefix>_class12.txt` rather than receiving an invented fallback category.

## Preliminary `class4` and detailed descriptions

`class4` is a diagnostic intermediate with these values:

- `new_junction`: at least one query junction does not match any same-locus
  reference junction within the offset;
- `all_matched>=_<ref>`: a splice chain matches and query genomic span is at
  least the reference span;
- `all_matched_<_<ref>`: a splice chain matches and query genomic span is
  shorter than the reference span;
- `new_combination`: no new junction exists, but no complete reference splice
  chain matches.

`class12` does not blindly copy the `all_matched` direction. It applies the
Figure 2I summed-exon 5% rule against the same selected reference.

The `_desc.txt` miss/extra fields retain detailed compatibility descriptions,
including 5′/3′ terminal missing/extra, intron retention, and internal exon
missing/extra.

## Optional end-shift diagnostics

`--end-shift-bp <N>` adds strand-aware diagnostic text for splice-equal records
to `_desc.txt`; `0` disables it and is the default.

- `5 end extension: <bp>` / `5 end truncation: <bp>`
- `3 end extension: <bp>` / `3 end truncation: <bp>`

This option never changes `class12`. Its base-pair cutoff controls only whether
the detailed text is printed and is independent of the mandatory Figure 2I 5%
classification threshold. The former Rust-only values
`5'end_extension`, `3'end_extension`, `5'end_truncation`, and
`3'end_truncation` are not Figure 2 categories and are no longer emitted.

## Important parameters

- `--offset-bp` defaults to `10` bp and controls fuzzy junction matching.
- `--fusion-fraction-read` defaults to `0.1`.
- `--fusion-fraction-ref` defaults to `0.1`.
- Fusion detection ignores the input gene field, reassigns genes by overlap,
  and requires two or more distinct genes after applying both fractions.
