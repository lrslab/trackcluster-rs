# m6A / modification calls isoform-level analysis plan

## 背景和目标

目标是在现有 `trackcluster-rs` 的 isoform clustering/counting 结果之上，增加对 modification calls 的 isoform 层面汇总和差异分析。输入 call 可能来自：

- Dorado/modBAM：BAM/SAM 中的 `MM`/`ML`/`MN` modified-base tags。
- m6Anet RNA002：`data.indiv_proba.csv` 和/或 `data.site_proba.csv`。

核心原则：modification calls 不应写入 BED `score` 或 `extra_fields`，而应作为 sidecar 表和最终 `read_to_isoform` 映射 join。这样可以兼容当前默认 `--name2-mode coverage`，也不会污染 TrackCluster 已有的 bigGenePred extra field 语义。

## 多代理代码审阅结论

1. 当前仓库是 BED12/bigGenePred-like 原生流程，没有 BAM/SAM/PAF/GTF/modification parser。
   - `src/io/mod.rs` 只导出 `bed` 和 `manifest`。
   - `src/io/bed.rs` 的 `read_bed12()` 是唯一结构化转录本/reads reader。
   - `Cargo.toml` 目前只有 `anyhow`、`clap`、`thiserror` 等轻依赖。

2. 核心数据模型是 `src/model/transcript.rs::Transcript`。
   - 包含 `chrom`、`strand`、`tx_start`、`tx_end`、`name`、`score`、`exons`、`extra_fields`。
   - 没有 BAM read sequence、CIGAR、query offset、base quality 或 modified-base tag 信息。

3. 当前 isoform 支撑读段关系的稳定接口是 `read_to_isoform`。
   - `clusterj` / `cluster` 返回 `ClusterResult { isoforms, read_to_isoform, unused }`。
   - CLI 和 `flow` 会写出 `*_read_to_isoform.tsv`。
   - unique counting 会额外写 `*_read_to_isoform.unique.tsv`，这是最终 count 使用的 read assignment。

4. 多样本流程会重写 read id。
   - manifest pooling 会把 read name 改成 `<sample>::<orig_read_id>`。
   - modification call 输入必须保留或可恢复这个 sample 前缀，否则不能可靠做 sample/group 层统计。

5. 现有 `count-multi` 是 expression/usage 专用，不适合直接塞 modification 差异统计。
   - 建议新增 `src/modification/` 模块，只复用 manifest、BED、read-to-isoform、sample tagging 等公共能力。

## 外部格式依据

- Dorado modified-base docs：Dorado modified-base output 写入 SAM/BAM `MM`/`ML` tags，`--modified-bases-threshold` 只影响 tags 中哪些概率被写出，不改变模型原始概率。
  - https://software-docs.nanoporetech.com/dorado/latest/basecaller/mods/
- SAM optional fields spec：`MM` 位置是按 as-sequenced read orientation 记录，需要结合 `ML` 概率、`MN` 长度、FLAG/CIGAR 投影后才能落到参考或 transcript 坐标。
  - https://samtools.github.io/hts-specs/SAMtags.pdf
- ONT BAM output spec：modified bases 存在时 `MN`、`ML`、`MM` 是相关字段。
  - https://nanoporetech.github.io/ont-output-specifications/latest/read_formats/bam/
- m6Anet docs：`data.indiv_proba.csv` 是 read-level probability；`data.site_proba.csv` 是 site-level probability/mod_ratio。m6Anet RNA002 模型有默认 read probability threshold。
  - https://m6anet.readthedocs.io/en/latest/

## 建议的用户可见能力

### 1. 标准化 modification call 输入

先定义一个内部 normalized TSV，作为 Dorado 和 m6Anet 的共同中间层：

```tsv
sample	read_id	source	mod_code	chrom	pos0	strand	transcript_id	transcript_pos0	kmer	probability	is_modified
S1	S1::read001	dorado	m6A	chr1	12345	+	NA	NA	DRACH	0.92	true
S1	S1::read002	m6anet	m6A	NA	NA	+	ENST0001	77	AGACT	0.18	false
```

字段语义：

- `sample`: manifest sample id；单样本可允许为空或自动填 `sample`。
- `read_id`: 必须和 `*_read_to_isoform.tsv` 中的 read id 完全一致；manifest 模式下应包含 `<sample>::`。
- `source`: `dorado`、`m6anet` 或 `normalized`。
- `mod_code`: 初期聚焦 `m6A`，但类型保留为字符串。
- `chrom` / `pos0` / `strand`: 0-based genomic site；Dorado 推荐填。
- `transcript_id` / `transcript_pos0`: transcript coordinate site；m6Anet 推荐填。内部统一为 0-based。
- `kmer`: 可选 motif/context。
- `probability`: 0.0 到 1.0。
- `is_modified`: 由 caller 原始判定或用户指定 threshold 得到；如果为空，可在聚合阶段按 `--prob-threshold` 计算。

### 2. Isoform-level aggregation

新增命令建议：

```bash
trackcluster mod-aggregate \
  --manifest samples.tsv \
  --isoform out/prefix_isoform.bed \
  --read-to-isoform out/prefix_read_to_isoform.unique.tsv \
  --mod-calls calls.normalized.tsv \
  --out out/prefix
```

输出：

```text
<prefix>.isoform_mod_sites.tsv
<prefix>.isoform_mod_burden.tsv
```

`isoform_mod_sites.tsv` 建议列：

```tsv
gene	isoform_id	sample	group	mod_code	site_key	chrom	pos0	strand	transcript_id	transcript_pos0	n_callable	n_modified	mod_ratio	mean_probability	kmer	source
```

`isoform_mod_burden.tsv` 建议列：

```tsv
gene	isoform_id	sample	group	mod_code	n_sites	n_callable	n_modified	mod_ratio	mean_probability
```

聚合规则：

- 用最终 `read_to_isoform` join modification call 到 isoform。
- 对 multi-mapped raw mapping，支持 `--assignment-mode fractional`；默认建议使用 `.unique.tsv`，避免同一个 read 同时贡献给多个 isoform。
- site key 优先用 genomic coordinate：`chrom:pos0:strand:mod_code`。
- 只有 transcript coordinate 时，用 `transcript_id:transcript_pos0:strand:mod_code`，并在可用参考注释时再投影到 genomic coordinate。
- `n_callable` 必须只统计该 site 上有明确概率/调用信息的 read；不能把未出现在 call 表中的 read 默认当 unmodified，除非输入格式明确提供低概率/negative calls。
- Dorado `MM` tag 受 threshold 影响，缺失 call 不一定是 unmodified；如果要估计 stoichiometry，建议要求 Dorado 导出 threshold 0 或使用能提供 per-base probabilities 的 normalized sidecar。

### 3. Isoform-level differential analysis

新增命令建议：

```bash
trackcluster mod-diff \
  --manifest samples.tsv \
  --aggregate out/prefix.isoform_mod_sites.tsv \
  --group-a control \
  --group-b treated \
  --out out/prefix
```

输出：

```text
<prefix>.isoform_mod_diff.tsv
<prefix>.isoform_mod_pairwise.tsv
```

`isoform_mod_diff.tsv` 建议列：

```tsv
gene	isoform_id	mod_code	site_key	group_a	group_b	n_callable_a	n_modified_a	ratio_a	n_callable_b	n_modified_b	ratio_b	delta_ratio	odds_ratio	p_value	q_value
```

第一阶段统计策略：

- 对每个 `(isoform_id, site_key, mod_code)` 做 group-level pooled 2x2 Fisher exact test。
- 用 Benjamini-Hochberg 计算 `q_value`。
- 增加过滤参数：`--min-callable-per-group`、`--min-samples-per-group`、`--min-delta-ratio`。
- 同时输出 per-sample aggregate，避免只看 pooled 结果造成 replicate 被掩盖。

第二阶段统计策略：

- 增加 replicate-aware beta-binomial 或 logistic mixed model 支持。
- 如果不想引入重依赖，可先输出设计矩阵，由 R/Python 下游完成复杂模型。

Isoform 间区别分析：

- 在同一 gene 内比较不同 isoform 的 shared genomic sites。
- 输出 `isoform_mod_pairwise.tsv`，列出 `(gene, site_key, isoform_a, isoform_b, ratio_a, ratio_b, delta_ratio, p/q)`。
- 对 isoform-specific exons/sites，单独标记为 `site_presence=isoform_specific`，不要和 shared site 混入同一统计检验。

## 代码分层计划

### Phase 0: 固化格式和边界

- 在 `docs/FORMATS.md` 增加 normalized modification calls TSV、aggregation output、diff output schema。
- 在 `docs/CLI.md` 增加 `mod-aggregate` / `mod-diff` 草案。
- 明确第一版必须支持 normalized TSV；Dorado/m6Anet parser 可分阶段接入。

### Phase 1: 核心类型和 TSV parser

新增模块：

```text
src/modification/mod.rs
src/modification/call.rs
src/modification/aggregate.rs
src/modification/diff.rs
src/modification/stats.rs
src/io/mod_calls.rs
src/cli/mod_aggregate.rs
src/cli/mod_diff.rs
```

主要类型：

```rust
struct ModificationCall {
    sample: Option<String>,
    read_id: String,
    source: ModificationSource,
    mod_code: String,
    site: ModificationSite,
    kmer: Option<String>,
    probability: Option<f64>,
    is_modified: Option<bool>,
}

enum ModificationSite {
    Genomic { chrom: String, pos0: Coord, strand: Strand },
    Transcript { transcript_id: String, pos0: u32, strand: Strand },
}
```

需要实现：

- normalized TSV reader/writer。
- probability threshold 应用。
- read id sample prefix 校验和可选补齐。
- join `ModificationCall` 与 `read_to_isoform`。
- 输出 per-sample isoform site aggregate。

### Phase 2: 坐标转换 helper

扩展 `Transcript`，但不要改变现有 BED 语义：

- `contains_genomic_pos(pos: Coord) -> bool`
- `genomic_to_spliced_offset(pos: Coord) -> Option<u32>`
- `spliced_offset_to_genomic(offset0: u32) -> Option<Coord>`
- strand-aware 规则：plus strand 从低坐标 exon 到高坐标累加；minus strand 从高坐标 exon 到低坐标累加。

用途：

- 将 m6Anet transcript coordinate 投影到 reference transcript genomic coordinate。
- 将 genomic site 转换成每个 isoform 的 spliced coordinate，方便判断 site 是否落在 isoform exon 内。
- 区分 shared genomic site、isoform-specific site、不可投影 site。

### Phase 3: m6Anet RNA002 adapter

新增 `src/io/m6anet.rs`：

- 支持 `data.indiv_proba.csv` read-level 输入，优先用于 isoform-level 分析。
- 支持 `data.site_proba.csv` site-level 输入，但标记为 site summary；如果没有 read id，不能做 read-to-isoform join，只能做 transcript/site 注释或参考 transcript summary。
- 支持 gzip CSV：可能需要 `flate2` 或让用户先解压；第一版可先只支持明文 CSV。
- 需要 read index 到 read id 的映射输入，例如 Nanopolish summary 或 m6Anet dataprep 关联表；没有映射就无法和 `read_to_isoform` join。
- RNA002 默认参数可以提供 `--source-preset m6anet-rna002`，设置推荐 `--prob-threshold`，但必须允许用户覆盖。

CLI 例子：

```bash
trackcluster mod-aggregate \
  --source m6anet \
  --m6anet-indiv data.indiv_proba.csv \
  --m6anet-read-index summary.txt \
  --manifest samples.tsv \
  --isoform pooled_isoform.bed \
  --read-to-isoform pooled_read_to_isoform.unique.tsv \
  --out pooled
```

### Phase 4: Dorado / modBAM adapter

两条路线：

1. 保守路线：先要求用户提供 normalized TSV。
   - 可以由外部工具从 Dorado BAM/CRAM 生成。
   - 这能先验证 isoform aggregation 和 differential analysis 逻辑。

2. 原生路线：新增 `src/io/modbam.rs`。
   - 依赖候选：`noodles` 或 `rust-htslib`。
   - 解析 BAM/SAM QNAME、FLAG、CIGAR、SEQ、`MM`、`ML`、`MN`。
   - 校验 `MN` 与 SEQ length，处理 soft/hard clipping、reverse-complement orientation。
   - 将 read-base modification calls 投影到 reference coordinate。
   - 如果 BAM 是 transcriptome alignment，再输出 transcript coordinate；如果是 genome alignment，则输出 genomic coordinate。

第一版建议优先做保守路线，避免 BAM tag/CIGAR 细节阻塞主功能。

### Phase 5: differential statistics

实现轻量统计：

- Fisher exact test。
- Benjamini-Hochberg FDR。
- coverage/sample filters。
- deterministic sorting，保证 golden test 稳定。

后续可选：

- beta-binomial overdispersion。
- per-gene or per-isoform covariates。
- 导出 R/Python design matrix。

### Phase 6: flow 集成

在 `src/flow/full.rs` 中，最终 `count_read_to_isoform` 确定后增加可选 modification aggregation：

- `flow --mod-calls calls.normalized.tsv`
- `flow --mod-source m6anet|dorado|normalized`
- `flow --mod-out-prefix <prefix>` 默认随 flow prefix。

注意：

- unique assignment 模式下应使用 `<prefix>_read_to_isoform.unique.tsv`。
- fractional 模式下用 raw `<prefix>_read_to_isoform.tsv` 并保留 fractional 权重。
- manifest 模式必须检查 modification call 的 sample/read id 与 pooled read id 一致。

## 测试计划

新增 fixtures：

```text
tests/fixtures/moddiff/
  isoforms.bed
  reads.bed
  read_to_isoform.unique.tsv
  samples_2x2.tsv
  calls.normalized.tsv
  m6anet_indiv.csv
  m6anet_read_index.tsv
```

新增 golden：

```text
tests/golden/moddiff/
  isoform_mod_sites.tsv
  isoform_mod_burden.tsv
  isoform_mod_diff.tsv
```

单元测试：

- normalized TSV parser：缺字段、空 probability、非法 probability、非法 coordinate。
- sample prefix：`S1::read` 与 manifest sample 匹配；未加前缀时按参数补齐。
- read-to-isoform join：missing read、missing isoform、duplicate calls。
- coordinate helper：plus/minus strand、多 exon、exon boundary、intron position。
- aggregation：same read multi-site、same site multi-read、threshold 后 `is_modified`。
- stats：Fisher exact、BH FDR、zero denominator、排序稳定。

集成测试：

- `trackcluster mod-aggregate` 从 normalized TSV 生成 aggregate golden。
- `trackcluster mod-diff` 从 aggregate 表生成 differential golden。
- 2 group x 2 replicate：至少覆盖显著差异、不显著、低覆盖过滤、方向相反。
- m6Anet adapter：read_index 映射到 read_id 后可 join。

## 关键风险和需要明确的决策

1. Dorado 缺失 call 的含义。
   - 如果 `--modified-bases-threshold` 大于 0，BAM 中缺失的低概率 call 不能简单当作 unmodified。
   - 计划中需要显式区分 `missing`、`callable_unmodified`、`modified`。

2. m6Anet site-level 文件缺少 read id。
   - `data.site_proba.csv` 可以做 transcript/site summary，但不能可靠做 read-to-isoform aggregation。
   - isoform-level read assignment 应优先要求 `data.indiv_proba.csv`。

3. transcript coordinate 的参考对象。
   - m6Anet `transcript_id` 指向参考 transcript，不一定等于 trackcluster 输出 isoform id。
   - 需要参考 transcript BED/GTF 来投影成 genomic site，再和发现的 isoform exon 结构关联。

4. read id 规范化。
   - Dorado QNAME、Nanopolish/m6Anet read index、TrackCluster pooled read id 必须建立一条可审计映射链。
   - manifest 模式下建议在 normalized calls 中强制使用 `<sample>::<orig_read_id>`。

5. 统计独立性。
   - 多 reads 来自同一样本，pooled Fisher 只能作为第一版探索性结果。
   - 有 biological replicates 时，应优先输出 per-sample rates，并为后续 replicate-aware 模型预留接口。

## 推荐实施顺序

1. 写 normalized TSV schema、`mod-aggregate`、`mod-diff`，只支持明文 TSV。
2. 加 coordinate helpers 和完整单元测试。
3. 用 synthetic normalized calls 建 fixtures/goldens。
4. 接 m6Anet `data.indiv_proba.csv` adapter。
5. 再评估 Dorado 原生 BAM parser；短期先接受外部 normalized TSV。
6. 最后把可选 aggregation 接入 `flow`。

