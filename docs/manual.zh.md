# rnaseq-standard-flow 中文用户手册

`rnaseq-standard-flow` 是 TAFFISH RNA-seq 标准流程的总入口。它把已经发布并固定版本的 RNA-seq 子流程串联起来，让用户可以从常见 bulk RNA-seq 输入出发，得到表达矩阵、差异表达结果、功能富集结果、图表、日志、版本记录和一份静态 HTML 项目报告。

本手册面向三类用户：

- 只想快速跑通一套标准 RNA-seq 分析的新用户。
- 需要知道每个参数和输入文件如何准备的生物信息用户。
- 需要审计结果、复跑命令和交付报告的项目维护者。

## 1. 流程定位

默认路线是轻量的 Salmon-first 路线：

```text
FASTQ + genome + annotation + metadata + GMT
-> rnaseq-index-flow
-> rnaseq-expression-flow
-> rnaseq-de-flow
-> rnaseq-enrichment-flow
-> rnaseq-report-flow
```

显式指定 `--route both` 时，流程同时运行 alignment/count/QC 分支：

```text
rnaseq-index-flow --genome-indexer hisat2
-> rnaseq-alignment-flow
-> rnaseq-alignment-qc-flow
-> rnaseq-count-flow
```

默认情况下，差异表达使用 Salmon/tximport 生成的 gene-level count matrix。使用 `--route both --de-source featurecounts` 时，差异表达改用 featureCounts 的 gene count matrix。

这个流程不是 Nextflow 或 Snakemake 式的大型 workflow engine。它是一个 TAFFISH 风格的 shell-native umbrella flow：每个核心分析步骤都通过固定版本的 `taf-rnaseq-*-flow` 子流程执行，并把命令、版本、日志和结果集中写入一个显式输出目录。

## 2. 最小快速开始

最常见的运行方式：

```sh
taf-rnaseq-standard-flow \
  --samples samples.tsv \
  --genome genome.fa \
  --annotation genes.gff3 \
  --metadata metadata.tsv \
  --design '~ condition' \
  --contrast condition:treated:control \
  --gene-sets gene_sets.gmt \
  --background background.tsv \
  --outdir rnaseq-standard-out
```

如果你还需要 HISAT2 BAM、BAM 层面 QC 和 featureCounts 计数，但差异分析仍使用 Salmon：

```sh
taf-rnaseq-standard-flow \
  --samples samples.tsv \
  --genome genome.fa \
  --annotation genes.gff3 \
  --metadata metadata.tsv \
  --design '~ condition' \
  --contrast condition:treated:control \
  --gene-sets gene_sets.gmt \
  --background background.tsv \
  --outdir rnaseq-standard-out \
  --route both
```

如果你希望差异分析使用 featureCounts count matrix：

```sh
taf-rnaseq-standard-flow \
  --samples samples.tsv \
  --genome genome.fa \
  --annotation genes.gff3 \
  --metadata metadata.tsv \
  --design '~ condition' \
  --contrast condition:treated:control \
  --gene-sets gene_sets.gmt \
  --background background.tsv \
  --outdir rnaseq-standard-out \
  --route both \
  --de-source featurecounts
```

## 3. 输入数据准备

### 3.1 FASTQ 样本表

`--samples` 指向一个 tab-separated 文件。最小 single-end 示例：

```text
sample_id	read1	condition
S1	data/S1.fastq.gz	control
S2	data/S2.fastq.gz	treated
```

最小 paired-end 示例：

```text
sample_id	read1	read2	condition
S1	data/S1_R1.fastq.gz	data/S1_R2.fastq.gz	control
S2	data/S2_R1.fastq.gz	data/S2_R2.fastq.gz	treated
```

规则：

- `sample_id` 必须唯一且非空。
- `read1` 必须存在。
- `read2` 存在时表示 paired-end。
- FASTQ 相对路径按 `samples.tsv` 所在目录解释，不按当前目录或输出目录解释。
- 默认 Salmon 路线可以处理 single-end 或 paired-end。
- `--route both` 要求所有样本都是同一种 layout。不能一部分 single-end、一部分 paired-end。

建议额外保留 `condition`、`batch`、`library_layout`、`strandedness` 等列，便于检查和后续报告，但本 flow 的正式设计信息以 `--metadata` 为准。

### 3.2 参考基因组 FASTA

`--genome` 指向参考基因组 FASTA，例如：

```text
>chrI
ACGT...
>chrII
ACGT...
```

注意：

- annotation 中的 sequence ID 必须和 FASTA header 的第一个 token 匹配。
- 如果 FASTA header 是 `>chrI some description`，annotation 应使用 `chrI`。
- 不建议在正式项目中混用不同数据库版本的 genome 和 annotation。

### 3.3 注释文件 GFF3/GTF

`--annotation` 支持 GFF3 或 GTF。它会由 `rnaseq-index-flow` 标准化，并用于：

- 提取 transcript FASTA。
- 生成 `tx2gene.tsv`。
- 生成 gene-level GTF/GFF3。
- 在 `--route both` 时供 featureCounts 和 RNA-seq QC 使用。

常见问题是 FASTA 与 annotation 的 seqid 不一致。例如 FASTA 使用 `chrI`，GFF3 使用 `I`。这种情况应在数据准备阶段统一。

### 3.4 metadata 表

`--metadata` 传给 `rnaseq-de-flow`。默认样本列名是 `sample`：

```text
sample	condition	batch
S1	control	B1
S2	control	B1
S3	treated	B2
S4	treated	B2
```

规则：

- 默认 `sample` 列必须和 count matrix 的样本名一致。
- 如果你的列名不是 `sample`，使用 `--sample-column`。
- `condition`、`batch` 等列可用于 DESeq2 design。
- 不要在 metadata 中混入 count matrix 没有的样本，反之亦然。

### 3.5 design 和 contrast

`--design` 是 DESeq2 design formula。例如：

```text
--design '~ condition'
```

如果有 batch：

```text
--design '~ batch + condition'
```

`--contrast` 使用 `FACTOR:NUMERATOR:DENOMINATOR`：

```text
--contrast condition:treated:control
```

含义是计算：

```text
treated / control
```

因此 log2 fold change 大于 0 表示 treated 组表达更高，小于 0 表示 control 组表达更高。

### 3.6 GMT gene set 文件

`--gene-sets` 指向离线 GMT 文件。格式：

```text
term_id	description	gene1	gene2	gene3
```

本 flow 不联网下载 GO、KEGG、Reactome 或其他数据库。gene set 文件必须由用户准备好，并且 gene ID 空间要和 DE 结果一致。

### 3.7 background gene universe

`--background` 是可选但推荐的背景基因列表。最常见格式：

```text
gene_id
YAL001C
YAL002W
```

提供 background 时，standard-flow 会在调用 enrichment 前，把 DE significant gene list 和 ranked gene table 过滤到这个 ID 空间，并记录：

```text
04_reports/enrichment_background_filter.tsv
```

这样可以避免少量 annotation-only 或 gene-set-only ID 破坏 enrichment 的背景一致性。

## 4. 常规运行模式

### 4.1 默认 Salmon-first 路线

适合：

- 常规 bulk RNA-seq 表达定量。
- 教学和轻量项目。
- 不需要 BAM 交付的项目。
- 希望资源使用较轻的快速分析。

默认路线输出 gene-level counts、TPM、DESeq2 结果、ORA/GSEA 结果和最终报告。

### 4.2 `--route both`

适合：

- 需要 sorted BAM。
- 需要 HISAT2 alignment summary。
- 需要 RSeQC / Qualimap / SAMtools 的 BAM 层面 QC。
- 需要保留传统 alignment-count route 证据。

这个模式会增加运行时间、磁盘和内存需求。默认差异分析仍使用 Salmon counts，除非同时设置 `--de-source featurecounts`。

### 4.3 `--route both --de-source featurecounts`

适合：

- 项目要求 classical genome alignment + featureCounts + DESeq2。
- 需要让差异分析完全基于 BAM/featureCounts。
- 需要把 Salmon 结果和 featureCounts 结果作为可比较的两条证据线。

注意：这不是“比 Salmon 更正确”的通用结论。Salmon-first 和 alignment-count 是不同建模路线，选择应取决于项目目标、样本质量、参考注释质量和交付需求。

## 5. 参数详解

### 5.1 必填参数

| 参数 | 含义 |
| --- | --- |
| `--samples` | FASTQ 样本表。 |
| `--genome` | 参考基因组 FASTA。 |
| `--annotation` | GFF3/GTF 注释。 |
| `--metadata` | DESeq2 metadata 表。 |
| `--design` | DESeq2 design formula。 |
| `--contrast` | 差异比较，格式为 `factor:numerator:denominator`。 |
| `--gene-sets` | 离线 GMT gene set 文件。 |
| `--outdir` | 专用输出目录。已有目录默认拒绝。 |

### 5.2 主要流程参数

| 参数 | 默认值 | 含义 |
| --- | --- | --- |
| `--threads` | `2` | 传给 index、expression 和相关子流程的线程数。 |
| `--route` | `salmon` | `salmon` 或 `both`。 |
| `--de-source` | `salmon` | `salmon` 或 `featurecounts`。 |
| `--project-name` | `RNA-seq project` | 最终 HTML 报告显示的项目名。 |
| `--force` | off | 覆盖已有 standard-flow 输出目录中的本流程产物。 |

### 5.3 expression 参数

| 参数 | 默认值 | 含义 |
| --- | --- | --- |
| `--library-type` | `A` | Salmon library type。`A` 表示自动推断。 |
| `--indexer` | `salmon` | reference 阶段建立 `salmon` 或 `both` 索引。r2 quant 仍使用 Salmon。 |
| `--kmer` | `31` | Salmon/Kallisto index k-mer 参数。 |
| `--trim` | off | 在 expression subflow 中启用 fastp trimming。 |
| `--skip-fastqc` | off | 跳过 expression subflow 中的 FastQC。 |
| `--min-assigned-frags` | `10` | Salmon quant 最少 assigned fragments 检查阈值。 |
| `--counts-from-abundance` | `no` | tximport gene count 处理方式：`no`、`scaledTPM`、`lengthScaledTPM`、`dtuScaledTPM`。 |

### 5.4 alignment/count/QC 参数

这些参数主要在 `--route both` 时有意义。

| 参数 | 默认值 | 含义 |
| --- | --- | --- |
| `--rna-strandness` | `none` | HISAT2 RNA strandness：`none`、`F`、`R`、`FR`、`RF`。 |
| `--alignment-min-mapq` | `0` | alignment 和 count 相关 MAPQ 过滤阈值。 |
| `--count-strand` | `0` | featureCounts strand mode：`0` unstranded、`1` stranded、`2` reversely stranded。 |
| `--count-feature-type` | `exon` | featureCounts 统计的 feature 类型。 |
| `--count-attribute` | `gene_id` | featureCounts 分组属性。 |
| `--count-min-assigned-reads` | `0` | count-flow assigned reads 最低检查阈值。 |
| `--qc-mapq` | `30` | RSeQC QC 使用的 MAPQ cutoff。 |
| `--infer-sample-size` | `200000` | RSeQC infer_experiment.py 抽样大小。 |
| `--java-mem-size` | `4G` | Qualimap Java memory 参数。 |
| `--sequencing-protocol` | `non-strand-specific` | Qualimap 协议：`non-strand-specific`、`strand-specific-forward`、`strand-specific-reverse`。 |

### 5.5 DE 参数

| 参数 | 默认值 | 含义 |
| --- | --- | --- |
| `--sample-column` | `sample` | metadata 中的样本列。 |
| `--gene-column` | `gene_id` | count matrix 的基因列。 |
| `--padj-cutoff` | `0.05` | 显著差异基因 adjusted P-value 阈值。 |
| `--lfc-cutoff` | `1` | 显著差异基因绝对 log2 fold change 阈值。 |
| `--fit-type` | `parametric` | DESeq2 dispersion fit：`parametric`、`local`、`mean`。 |
| `--lfc-shrink` | `none` | LFC shrinkage：`none`、`ashr`、`apeglm`。 |
| `--coef` | empty | `--lfc-shrink apeglm` 时需要指定 DESeq2 coefficient。 |
| `--min-count` | `1` | 低表达基因过滤的最低 count。 |
| `--min-samples` | `2` | 满足最低 count 的最少样本数。 |
| `--top-var` | `500` | PCA/样本结构图中选择的 top variable genes 数。 |
| `--top-heatmap` | `50` | heatmap 展示的 top genes 数。 |

### 5.6 enrichment 参数

| 参数 | 默认值 | 含义 |
| --- | --- | --- |
| `--enrichment-min-size` | `2` | GMT gene set 最小大小。 |
| `--enrichment-max-size` | `500` | GMT gene set 最大大小。 |
| `--enrichment-pvalue-cutoff` | `1` | enrichment P-value cutoff。默认保留完整结果供报告和筛选。 |
| `--enrichment-padj-method` | `BH` | 多重检验校正方法。 |
| `--enrichment-top-n` | `20` | 富集图展示的 top gene sets 数。 |
| `--enrichment-seed` | `1` | GSEA 等步骤使用的随机种子。 |

## 6. 输出目录和关键文件

标准输出结构：

```text
rnaseq-standard-out/
  00_inputs/
  01_logs/
  02_intermediate/
  03_results/
    reference/
    expression/
    alignment/
    alignment_qc/
    count/
    de/
    enrichment/
    report/
    plots/
      png/
      pdf/
  04_reports/
    rnaseq_report.html
    report_interpretation.html
    plot_files.tsv
    commands.sh
    versions.tsv
    methods.txt
    flow_summary.tsv
    subflows.tsv
  run.manifest.json
```

重点文件：

| 文件 | 用途 |
| --- | --- |
| `04_reports/rnaseq_report.html` | 最终项目报告入口。 |
| `04_reports/report_interpretation.html` | 配套 RNA-seq 解读指南，包含更长的生物学和技术解释。 |
| `03_results/report/` | report-flow 的完整收集结果。 |
| `03_results/plots/png/` | 顶层 PNG 图集。 |
| `03_results/plots/pdf/` | 顶层 PDF 图集。 |
| `04_reports/plot_files.tsv` | 顶层图集索引。 |
| `04_reports/commands.sh` | 主流程与上游子流程命令记录。 |
| `04_reports/versions.tsv` | 版本记录。 |
| `04_reports/methods.txt` | methods 文本。 |
| `04_reports/subflows.tsv` | 每个子流程状态。 |
| `run.manifest.json` | 机器可读运行 manifest。 |

## 7. 最终 HTML 报告解读

最终报告是静态 HTML，可直接在浏览器打开。它包含：

- 一键中文/英文切换。
- 左侧目录和当前章节高亮。
- 项目概览和关键指标。
- RNA-seq 技术路线图。
- reference、expression、alignment/count/QC、DE、enrichment、deliverables、provenance 分章节展示。
- PNG/PDF 图卡。
- MultiQC、FastQC、Qualimap 等 HTML 子报告链接。
- 工具和 TAFFISH source links。

常见图表含义：

| 图 | 含义 |
| --- | --- |
| PCA plot | 样本是否按条件、批次或异常样本结构聚类。 |
| Sample correlation heatmap | 样本间表达相关性。 |
| Expression distribution | 表达分布是否存在极端偏移。 |
| MA plot | 平均表达量与 log2 fold change 的关系。 |
| Volcano plot | 效应量和显著性的联合展示。 |
| DEG counts barplot | 上调/下调差异基因数量。 |
| Heatmap | top 差异基因在样本间的表达模式。 |
| Top genes expression | 代表性差异基因的表达情况。 |
| ORA dotplot/barplot | 显著基因列表中过度代表的 gene sets。 |
| GSEA NES plot | 排序基因列表中的方向性 gene set 信号。 |
| GSEA enrichment curves | 代表性 gene sets 在排序列表中的 running enrichment。 |

报告帮助用户理解结果结构，但不替代项目特异的生物学解释。最终结论应结合实验设计、样本质量、批次、验证实验和领域知识。

## 8. 使用 yeast 测试数据

可以用 `taf-rnaseq-yeast-get-data` 把 yeast SNF2 示例数据准备到你选择的本地目录。完整数据准备会下载公开 ENA/SGD/GO 资源，建议先确认网络和磁盘空间。

```sh
taf-rnaseq-yeast-get-data \
  --outdir yeast-snf2-data-v1 \
  --stage all \
  --resume true
```

数据准备完成后，把 `<outdir>/03_results` 当作数据根目录使用。下面的示例假设当前目录下已经有：

```text
yeast-snf2-data-v1/03_results/
```

可以先设置一个变量，减少命令长度：

```sh
DATA="$PWD/yeast-snf2-data-v1/03_results"
```

然后运行 standard-flow：

```sh
taf-rnaseq-standard-flow \
  --samples "$DATA/yeast-snf2-fastq-mini-v1/samples.tsv" \
  --genome "$DATA/yeast-reference-sgd-r64.4.1-v1/reference/genome/yeast_s288c_reference_genome_R64-4-1.fa" \
  --annotation "$DATA/yeast-reference-sgd-r64.4.1-v1/reference/annotation/yeast_s288c_gene_annotation_R64-4-1.gff3" \
  --metadata "$DATA/yeast-snf2-fastq-mini-v1/metadata.tsv" \
  --design '~ condition' \
  --contrast condition:snf2_KO:WT \
  --gene-sets "$DATA/yeast-sgd-go-gene-sets-r64.4.1-v1/gene_sets/sgd_go_bp.gmt" \
  --background "$DATA/yeast-sgd-go-gene-sets-r64.4.1-v1/background/yeast_background_genes.tsv" \
  --outdir yeast-standard-out \
  --route both \
  --de-source featurecounts
```

如果只想快速确认下载计划，可以先运行 `--stage plan`；如果只需要 reference 或 gene sets，也可以分别运行 `--stage reference` 或 `--stage genesets`。实际项目中，`metadata.tsv` 必须与所选样本完全对应；如果你手动筛选了样本子集，也要同步筛选 metadata。

## 9. 推荐参数选择

常规项目建议：

- 默认先用 Salmon route 跑通，确认 FASTQ、reference、metadata、GMT 都能对齐。
- 如果项目需要 BAM 或 alignment QC，再使用 `--route both`。
- 如果交付要求 classical counts，用 `--route both --de-source featurecounts`。
- 有可靠 gene universe 时提供 `--background`。
- 不确定 Salmon library type 时保留 `--library-type A`。
- 已经做过外部 reads cleaning 时，可以不用 `--trim`。
- 正式项目不建议使用过小的 `--top-var` 或 `--top-heatmap`，除非只是快速测试。

## 10. 常见问题

### annotation 和 genome seqid 不匹配

症状：index-flow 在 reference 检查或 transcript 提取阶段失败。

处理：统一 FASTA header 第一 token 与 GFF3/GTF seqid。

### metadata 样本不匹配

症状：DE flow 报告 metadata 和 count matrix 样本不一致。

处理：检查 `sample_id`、metadata 的 `sample` 列，以及是否使用了 `--sample-column`。

### GMT / background / DE gene ID 不一致

症状：enrichment 结果很少、为空，或 background filtering 后 gene 数很低。

处理：确认 count matrix、annotation、GMT、background 使用同一 gene ID 体系。

### mixed single-end 和 paired-end

症状：`--route both` 直接拒绝运行。

处理：拆成两个独立项目运行，或在数据准备阶段统一 layout。

### 输出目录已存在

症状：流程拒绝已有 `--outdir`。

处理：换新目录，或确认可以覆盖后使用 `--force`。

### report 链接打不开

优先从 `04_reports/rnaseq_report.html` 打开最终报告，不要只复制单个 HTML 文件离开整个输出目录。HTML 子报告和图像链接依赖同一个 output tree。

## 11. 边界

`rnaseq-standard-flow` 不做这些事情：

- 不联网下载 reference、annotation、gene sets 或数据库。
- 不自动选择生物学 design。
- 不自动判断 strandedness。
- 不自动删除失败样本。
- 不替代完整生产级 workflow engine。
- r2 可以构建 Kallisto index，但表达定量仍走 Salmon。

如果项目需要复杂 batch correction、多因素交互设计、splicing analysis、fusion detection、allele-specific expression、single-cell RNA-seq 或临床级报告，应在本流程输出的基础上设计额外分析。
