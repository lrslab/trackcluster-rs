# Figure 1. TrackCluster-RS connects long-read isoform discovery with sample-level quantification and optional RNA modification analysis

**a, Workflow.** Genome-aligned long reads and reference transcript annotations are converted to compatible tracks. Reads from one or more samples are pooled, partitioned by gene and clustered to construct a shared catalog of known and novel isoforms. Catalog structures support transcript interpretation; the final unique read-to-isoform mapping is shared by expression quantification and optional modification aggregation. Alignment and modification calling are upstream of TrackCluster-RS.

**b, Junction and terminal evidence.** The default junction-chain mode corrects low-support splice sites using read and reference support and collapses compatible truncated reads while retaining their membership. Independently supported alternative transcript ends can be retained. Optional spliced-leader (SL) evidence protects sufficiently supported alternative 5′ ends when enabled; the default is no SL protection. Supported alternative 3′ endpoints are protected within the same splice chain. Both forms of terminal support are determined from the corrected locus before batching. The cartoons depict related decisions, rather than separate passes in which terminal evidence is evaluated only after collapsing. Two-pass exon/intron overlap clustering is available as an alternative mode.

**c, Transcript interpretation.** Query structures are compared with a reference to identify terminal exon changes, UTR changes, exon inclusion/skipping, intron retention, alternative splice sites and gene fusions. These examples summarize the legacy classification contract of 11 event labels plus the annotated-reference label. UTR changes require a matched splice chain and an absolute summed-exon-length difference of at least 5% of the selected reference.

**d, Multi-sample quantification.** Sample identities are retained after pooled discovery. Unique assignment produces isoform-by-sample counts, within-gene isoform usage and optional group summaries. Each assigned molecule contributes to one isoform; when the flow subsamples genes, abundance outputs incorporate the corresponding scaling. Matrix shades illustrate the output layout and are not observed expression measurements.

**e, Optional RNA modification analysis.** External Dorado or m6Anet calls are normalized to genomic sites and joined to the same final unique assignments. Paired circles depict example per-read observations grouped by assigned isoform, not a binary modification state of an entire isoform. Modified, unmodified and unknown observations are distinguished from structural absence of a site. Where the denominator is defined, the modified fraction is the number of modified molecules divided by the number of callable molecules. Coverage and eligibility QC, Wilson intervals and descriptive effect-size contrasts are available. Caller/model/chemistry strata remain separate; unknown observations are not counted as unmodified and contrast p/q values are not reported as inferential results.

All tracks are schematic plus-strand examples; repeated reads do not specify numerical support thresholds. Implementation rules are strand-aware. I1, I2 and I3 are illustrative isoform identifiers. Orange 5′ ticks and 3′ endpoint dots in panel b highlight terminal evidence. No panel contains benchmark or experimental results.

## 中文说明

本图以“联合发现 - 共享转录本目录 - 逐样本定量”为主线，沿用旧 TrackCluster 的外显子块和读段轨道表达。a 展示完整流程，并明确表达定量与修饰聚合复用同一份唯一读段分配。b 展示剪接位点校正、兼容截短读段合并，以及有独立支持的 5′/3′ 末端保留；SL 证据为可选，默认不启用。c 保留旧图的结构事件分类。d 展示多样本计数、基因内异构体使用比例及组汇总。e 展示外部修饰调用结果如何关联到异构体，区分未修饰、未知和结构上不存在的位点。

配色中的青色、蓝色和赭色在 a、d、e 中分别对应示例异构体 I1、I2、I3；c 中黑色为参考、蓝色为查询。所有图形和矩阵深浅均为示意，没有使用虚构的实验结果。图采用英文标注，便于作为论文方法概览图继续修改。

## Visual references

- Chen Y et al. *Context-aware transcript quantification from long-read RNA-seq data with Bambu.* Nature Methods 20, 1187-1195 (2023). [Figure 1](https://www.nature.com/articles/s41592-023-01908-w/figures/1), [article](https://doi.org/10.1038/s41592-023-01908-w). Used as a reference for a clear discovery-assignment-quantification narrative; the figure here does not depict Bambu's NDR classifier or EM quantification.
- Li R et al. *Direct full-length RNA sequencing reveals unexpected transcriptome complexity during Caenorhabditis elegans development.* Genome Research 30, 287-298 (2020). [Figure 2 and article](https://pmc.ncbi.nlm.nih.gov/articles/PMC7050527/#GR251512LIF2), [DOI](https://doi.org/10.1101/gr.251512.119). Used as a reference for read/exon cartoons and transcript-event interpretation. All artwork was drawn anew.

## Implementation basis and exports

Checked against TrackCluster-RS version 0.3.1, commit `980dd56`, and the local documentation at creation (2026-09-25):

- `README.md` and `docs/PIPELINE.md`: pooled discovery, unique assignment, expression and optional modification flow.
- `docs/behavior/cluster.md`: junction correction, terminal evidence and alternative overlap mode.
- `docs/behavior/desc.md`: the 11-event plus reference classification contract.
- `docs/FORMATS.md` and `src/modification/aggregate.rs`: callable denominators, site states, QC and descriptive contrasts.

Files:

- `docs/figures/trackcluster_rs_fig1.svg`: editable vector text and shapes.
- `output/pdf/trackcluster_rs_fig1.pdf`: vector PDF, 183 × 204.8 mm, embedded Arial fonts.
- `docs/figures/trackcluster_rs_fig1_600dpi.png`: 600 dpi raster export from the PDF.
- `docs/figures/trackcluster_rs_fig1_preview.png`: compact preview.
- `docs/figures/build_overview.py`: reproducible figure source; requires ReportLab and Arial regular/bold. Font paths can be changed for other operating systems.

Rebuild from the repository root:

```bash
python3 docs/figures/build_overview.py
pdftoppm -r 600 -png -singlefile output/pdf/trackcluster_rs_fig1.pdf docs/figures/trackcluster_rs_fig1_600dpi
pdftoppm -scale-to 2100 -png -singlefile output/pdf/trackcluster_rs_fig1.pdf docs/figures/trackcluster_rs_fig1_preview
```

The figure is intentionally an overview of implemented functionality, not a full list of command-line options. Default unique assignment is illustrated; alternative fractional expression counting is omitted. Modification statistics describe read-derived quantities and do not treat reads or technical partitions as biological replicates.
