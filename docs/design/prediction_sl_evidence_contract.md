# 新预测模块的 SL evidence 输入契约

状态：已修正无 SL 证据时的 BAM 转换默认分数；新预测模块的类型约束仍为设计，
当前仓库尚无独立预测模块，也没有实现下述 evidence 接口。

## 已核实的现状

- **没有 SL 证据的 read，其 BED score 默认为 `0`。**
  [`record_to_transcript`](../../src/io/bam.rs) 不提取 SL 证据，因此统一输出
  `Transcript.score = 0`，不随 MAPQ 大小变化。
- BAM `mapping_quality()` 仅用于 `--score` / `--min-mapq` 的比对质量过滤，
  缺失 MAPQ 按 `0` 过滤；原始 MAPQ 保留在 BAM 中，不占用 BED score。
- 转换只从 CIGAR 提取参考坐标与 exon blocks。softclip 不推进参考坐标，
  其长度、序列、质量值及原始 CIGAR 不写入
  [`Transcript`](../../src/model/transcript.rs)。比对跨度的端点仍保留，
  但这些坐标不是完整的原始 read 端点或 SL 命中证据。
- 转换器没有计算 SW score，也没有导入 SLRanger 的 `SL_score`。
- [`clusterj::is_sl_supported_read`](../../src/cluster/clusterj.rs) 在启用
  `--sw-score` 时直接比较 `tx.score > cutoff`；
  [`cluster_overlap::should_drop_read`](../../src/cluster/cluster_overlap.rs)
  的第二 pass 使用 `tx.score >= cutoff`。这些旧入口没有 score 来源类型校验。

因此，`bam2bigg → BED → --sw-score 11` 只会得到 score 为 `0`、在该阈值下
没有 SL 保护的 reads，不会从高 MAPQ 产生 SL 支持。该参数只解释已有 BED
score，不提取或验证 SL 序列证据。

当前 `clusterj`、`clusterj_batch` 和 `flow` 的 `--sw-score` 默认值为 `-1`。
直接运行 overlap `cluster` 的默认值为 `11`；当前转换输出的 `0` 分不达到
该阈值。这里的修复是转换器不再把 MAPQ 写入 BED score。

旧转换器已经生成的 MAPQ-scored BED 不会自动被修正，需要从 BAM 重新转换；
若继续使用这些旧文件，可用 `--sw-score -1` 禁用其 score 的 SL 解释。

## 独立、带类型的字段

新预测输入必须在 BED 几何模型之外保留独立的 `sl_evidence`，并将 MAPQ、
SLRanger `SL_score`、SW score 分别保存。建议的接口形状如下；这是设计示意，
不是已存在的 Rust API：

```text
PredictionRead {
    transcript: Transcript,
    mapq: Option<Mapq>,
    alignment_evidence: Option<AlignmentEvidence>,
    sl_evidence: SlEvidence,
}

SlEvidence {
    slranger: Option<SlRangerEvidence>,  // sl_score: SlRangerScore
    sw: Option<SwEvidence>,            // score: SmithWatermanScore
}
```

约束：

1. `Mapq`、`SlRangerScore`、`SmithWatermanScore` 必须是不同的类型，不能只是
   同一数值类型的别名，也不能提供跨类型的隐式转换或通用 score 回退逻辑。
   SLRanger 和 SW 证据可以同时存在，不能由一个字段相互覆盖。
2. 保留各来源原始值及其精度，不套用 BED 的 `0..=1000` 范围，不把
   `SL_score` 当作 SW score，不假定两者共享阈值 `11`。
   数值类型、范围及阈值规则须在导入器对接实际来源格式后确定。
3. 每份 SL 证据记录来源、工具/版本、原始记录身份及评价所需的上下文。
   SW 证据还须标明目标 SL 序列/模型、比对方向、命中区间和评分参数；
   一般序列比对的 SW 高分本身不能自动变成 SL 支持。
4. `sl_evidence` 缺失时是“未评估/无可用证据”，对应无 SL 的 BED 默认分数为
   `0`；这个默认值不表示已测得的 SW 或 SLRanger 分数。不能由高 MAPQ、BED score
   或 read 数量补足，也不能混同于“已评估但未获得支持”。导入失败应单独报告。
5. 分数阈值与证据类型绑定。SL 支持判定只消费经过验证的 SL 证据；MAPQ
   只参与比对质量筛选。禁止在新预测路径中通过重写 `Transcript.score` 或
   调用旧 score 判定来绕过这一约束。

## softclip、端点与持久化

需要端点序列证据的预测入口必须读取原始 BAM/序列或具有等价信息的证据文件，
不能从转换后的 BED 反推已经丢失的数据。

- `AlignmentEvidence` 应保留原始 CIGAR、可用的 query 序列/质量信息、左右
  softclip 区间以及原始比对端点，并明确坐标体系和序列方向。生物学 5′/3′
  必须按链及所用序列方向解释，不能一律等同于 CIGAR 左/右端。
- 明确区分“已观察到长度为零的 softclip”“缺失序列/质量信息”和 hardclip
  导致的不可观测部分；不能将这些状态统一写成空串或零。
- 每条证据关联 sample、read 及 alignment instance。重复 query name、次级或
  补充比对不能静默关联到另一条记录，也不能重复增加同一原始 read 的支持。
- 聚类或更换代表时保留证据来源及原始端点；吸收普通 reads 不能增加独立 SL
  测量数，也不能把证据位置移动到新代表端点。
- 若证据需要跨进程/断点续跑，使用带 schema 版本的独立文件或明确版本化的
  扩展格式；普通 BED12+8 不是这些证据的无损载体。记录来源及参数，并把证据
  文件和参数纳入结果的缓存身份。
- 普通 BED 输入默认没有 typed SL evidence。兼容既有 SL/SW-scored BED 时，
  必须由显式声明来源、分数语义和评价规则的适配器导入；仅设置旧 `--sw-score`
  不构成新模块的证据声明。来源不明的 BED score 也不能直接填入 `mapq`。

## 转换回归与后续验收

[`integration_aux_cli`](../../tests/integration_aux_cli.rs) 覆盖以下转换约束：

- MAPQ `60` 的 read 输出 BED score `0`；MAPQ `20` 仍被默认过滤阈值 `30`
  排除。通过 `--score 0` 或 `--min-mapq 0` 保留低 MAPQ read 后，BED score
  仍为 `0`。
- 正负链的高 MAPQ（`60`、`254`）、带 softclip reads 经过实际
  `bam2bigg → clusterj/cluster` 路径，在默认阈值及显式 `--sw-score 11` 下
  都不凭 MAPQ 获得替代 5′ 结构保护；结构、read 映射和 unused 均检查。
- 流式 BED 输出与内存 `read_bam` 的结果一致。

这些回归验证无 SL 的默认分数和过滤行为，不代表已提取端点序列证据，也不能
替代新预测模块的 typed evidence 验收。

新模块实现时必须加入以下验收：

- 高 MAPQ、没有 SL 证据的 reads 不产生 SL 支持；只有 BED 的输入也不产生。
- MAPQ、SLRanger `SL_score`、SW score 同时存在时，三个值及其来源分别往返
  保存；改变其中一个不能覆盖或重新解释另一个。
- SLRanger 与 SW 分别使用各自评价规则，并覆盖阈值边界；MAPQ 类型不能传给
  SL/SW 阈值接口。
- 缺失、未评估、已评估未支持和无效证据可区分；禁止缺失时回退到 BED score。
- 正负链、两端 softclip、hardclip、缺失序列及重复 alignment 身份都有覆盖。
- 合并、分批与持久化保留原始证据和端点；普通成员不被升级为 SL 测量支持。

本次修复不改变聚类阈值及其比较规则；实现新预测模块前，不能宣称现有 `Transcript`
已具备 typed SL evidence 或完整的端点序列证据。
