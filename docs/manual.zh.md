# rnaseq-standard-flow 中文用户手册

`rnaseq-standard-flow` 是 TAFFISH RNA-seq 标准流程的总入口。它把已经发布并固定版本的 RNA-seq 子流程串联起来，让用户可以从常见 bulk RNA-seq 输入出发，得到表达矩阵、差异表达结果、功能富集结果、图表、日志、版本记录和一份静态 HTML 项目报告。

本手册面向三类用户：

- 只想快速跑通一套标准 RNA-seq 分析的新用户。
- 需要知道每个参数和输入文件如何准备的生物信息用户。
- 需要审计结果、复跑命令和交付报告的项目维护者。

## 1. 流程定位

`0.3.0-r1` 默认仍是轻量的 reference / Salmon-first 路线，并兼容此前的
reference 用法：

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

显式指定 `--mode denovo` 时，流程进入无参路线：

```text
FASTQ + metadata + local protein FASTA + protein-to-GO map
-> rnaseq-denovo-assembly-flow
-> rnaseq-denovo-expression-flow
-> rnaseq-denovo-annotation-flow
-> rnaseq-de-flow
-> rnaseq-enrichment-flow
-> rnaseq-report-flow
```

无参路线不会因为用户漏传 `--genome` 或 `--annotation` 自动启动，必须显式设置
`--mode denovo`。这可以避免把原本应该做有参分析的数据误跑成无参分析。

如果项目主要使用无参路线，也建议阅读专门的无参分析手册：

- `docs/manual-denovo.zh.md`

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

如果没有可靠参考基因组和注释，并且已经准备好本地蛋白 FASTA 与 protein-to-GO
映射，可以使用无参模式：

```sh
taf-rnaseq-standard-flow \
  --mode denovo \
  --samples samples.tsv \
  --metadata metadata.tsv \
  --design '~ condition' \
  --contrast condition:treated:control \
  --protein-db proteins.faa \
  --go-map protein_go_map.tsv \
  --outdir rnaseq-denovo-standard-out
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

### 3.2 参考基因组 FASTA（reference 模式）

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

### 3.3 注释文件 GFF3/GTF（reference 模式）

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

### 3.6 GMT gene set 文件（reference 模式）

`--gene-sets` 指向离线 GMT 文件。格式：

```text
term_id	description	gene1	gene2	gene3
```

本 flow 不联网下载 GO、KEGG、Reactome 或其他数据库。gene set 文件必须由用户准备好，并且 gene ID 空间要和 DE 结果一致。

### 3.7 background gene universe（reference 模式）

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

### 3.8 无参模式的蛋白数据库和 GO 映射

`--mode denovo` 不需要 `--genome`、`--annotation`、`--gene-sets` 或
`--background`。它需要：

```text
--protein-db proteins.faa
--go-map protein_go_map.tsv
```

`--protein-db` 是本地蛋白 FASTA。`rnaseq-denovo-annotation-flow` 会用它建立
本地 DIAMOND 数据库，并把组装转录本预测出的蛋白与该数据库比对。

`--go-map` 是 protein subject ID 到 GO term 的映射表，推荐格式：

```text
subject_id	go_id	go_name	namespace
PROT1	GO:0006412	translation	biological_process
PROT2	GO:0005737	cytoplasm	cellular_component
```

无参模式会根据 DIAMOND 命中的 subject ID 和这张 GO map 生成：

```text
03_results/denovo_annotation/03_results/gene_sets/denovo_go.gmt
03_results/denovo_annotation/03_results/gene_sets/denovo_background.tsv
```

然后再把它们传给 enrichment 子流程。这里的 ID 是 assembled transcript ID，不是已知参考基因 ID；同源注释提供的是功能证据，不等同于人工审定注释。

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

### 4.4 `--mode denovo`

适合：

- 研究对象没有可靠参考基因组。
- 有参考基因组但注释质量很差，项目更关注转录本发现和表达模式。
- 需要为非模式物种先构建 transcriptome-level 表达矩阵。

无参模式会先组装转录本，然后对 assembled transcripts 做 Salmon 定量。差异表达使用 transcript-level count matrix，并固定 `--gene-column transcript_id`。这意味着结果中的 feature 不是 reference gene，而是组装转录本。后续富集依赖 `--protein-db` 和 `--go-map` 的同源注释质量。

无参模式不能使用 `--route both` 或 `--de-source featurecounts`，因为 BAM alignment、featureCounts 和 alignment QC 都需要参考基因组与注释。

## 5. 参数详解

### 5.1 必填参数

| 参数 | 含义 |
| --- | --- |
| `--samples` | FASTQ 样本表。 |
| `--mode` | `reference` 或 `denovo`。默认 `reference`，保持此前 reference 路线兼容。 |
| `--genome` | reference 模式必填；参考基因组 FASTA。 |
| `--annotation` | reference 模式必填；GFF3/GTF 注释。 |
| `--metadata` | DESeq2 metadata 表。 |
| `--design` | DESeq2 design formula。 |
| `--contrast` | 差异比较，格式为 `factor:numerator:denominator`。 |
| `--gene-sets` | reference 模式必填；离线 GMT gene set 文件。 |
| `--protein-db` | denovo 模式必填；本地蛋白 FASTA。 |
| `--go-map` | denovo 模式必填；protein subject ID 到 GO term 的 TSV。 |
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
| `--indexer` | `salmon` | reference 阶段建立 `salmon` 或 `both` 索引；reference 表达定量仍使用 Salmon。 |
| `--kmer` | `31` | Salmon/Kallisto index k-mer 参数。 |
| `--trim` | off | 在 expression subflow 中启用 fastp trimming。 |
| `--skip-fastqc` | off | 跳过 expression subflow 中的 FastQC。 |
| `--min-assigned-frags` | `10` | Salmon quant 最少 assigned fragments 检查阈值。 |
| `--counts-from-abundance` | `no` | tximport gene count 处理方式：`no`、`scaledTPM`、`lengthScaledTPM`、`dtuScaledTPM`。 |

### 5.4 denovo 参数

这些参数只在 `--mode denovo` 时有意义。

| 参数 | 默认值 | 含义 |
| --- | --- | --- |
| `--assembler` | `trinity` | 组装器：`trinity` 或 `rnaspades`。 |
| `--max-memory` | `4G` | 组装器内存限制，例如 `4G`、`32G`。 |
| `--min-contig-len` | `200` | 保留 assembled transcript 的最小长度。 |
| `--ss-lib-type` | `none` | strand-specific library type：`none`、`F`、`R`、`FR`、`RF`。 |
| `--no-normalize` | off | 通过 assembly-flow 传递 Trinity `--no_normalize_reads`。 |
| `--denovo-min-orf-aa` | `50` | TransDecoder 长 ORF 最小氨基酸长度。 |
| `--denovo-evalue` | `1e-5` | DIAMOND blastp e-value cutoff。 |
| `--denovo-max-target-seqs` | `1` | 每个预测蛋白保留的 DIAMOND subject hit 数。 |

正式项目中，`--max-memory` 应按照测序深度和物种转录组复杂度设置；小测试可以使用较小值。`--no-normalize` 主要用于已经下采样的小数据或测试数据，不是所有真实项目的默认最佳选择。

### 5.5 alignment/count/QC 参数

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

### 5.6 DE 参数

| 参数 | 默认值 | 含义 |
| --- | --- | --- |
| `--sample-column` | `sample` | metadata 中的样本列。 |
| `--gene-column` | `gene_id` | reference 模式 count matrix 的基因列；denovo 模式固定为 `transcript_id`。 |
| `--padj-cutoff` | `0.05` | 显著差异基因 adjusted P-value 阈值。 |
| `--lfc-cutoff` | `1` | 显著差异基因绝对 log2 fold change 阈值。 |
| `--fit-type` | `parametric` | DESeq2 dispersion fit：`parametric`、`local`、`mean`。 |
| `--lfc-shrink` | `none` | LFC shrinkage：`none`、`ashr`、`apeglm`。 |
| `--coef` | empty | `--lfc-shrink apeglm` 时需要指定 DESeq2 coefficient。 |
| `--min-count` | `1` | 低表达基因过滤的最低 count。 |
| `--min-samples` | `2` | 满足最低 count 的最少样本数。 |
| `--top-var` | `500` | PCA/样本结构图中选择的 top variable genes 数。 |
| `--top-heatmap` | `50` | heatmap 展示的 top genes 数。 |

### 5.7 enrichment 参数

| 参数 | 默认值 | 含义 |
| --- | --- | --- |
| `--enrichment-min-size` | `2` | GMT gene set 最小大小。 |
| `--enrichment-max-size` | `500` | GMT gene set 最大大小。 |
| `--enrichment-pvalue-cutoff` | `1` | enrichment P-value cutoff。默认保留完整结果供报告和筛选。 |
| `--enrichment-padj-method` | `BH` | 多重检验校正方法。 |
| `--enrichment-top-n` | `20` | 富集图展示的 top gene sets 数。 |
| `--enrichment-seed` | `1` | GSEA 等步骤使用的随机种子。 |

### 5.8 高级内部工具参数桥接

正常项目不需要使用这一节。`rnaseq-standard-flow 0.3.0-r1` 已经把常规运行需要的
关键参数作为顶层参数暴露出来。对于少数特殊项目，如果某个子流程内部工具需要额外参数，
standard-flow 会用带命名空间的结构域参数桥接到子流程内部的 `@step:`。

命名模式是：

```text
@<standard-step>-<child-tool-step>: ... @:
```

例如，`@denovo-assembly-trinity-assembly-step:` 会传给
`rnaseq-denovo-assembly-flow` 内部的 `@trinity-assembly-step:`。所有这些参数默认都是
空的，不会改变旧命令行为。完整槽位清单见 `docs/help.md` 和 README。常见分组如下：

| 分组前缀 | 控制的子流程内部步骤 |
| --- | --- |
| `@index-...` | `rnaseq-index-flow` 内部 AGAT/gffread/indexer 步骤 |
| `@expression-...` | `rnaseq-expression-flow` 内部 QC/修剪/Salmon/tximport/MultiQC 步骤 |
| `@alignment-...` | `rnaseq-alignment-flow` 内部 fastp/HISAT2/samtools/MultiQC 步骤 |
| `@alignment-qc-...` | `rnaseq-alignment-qc-flow` 内部 samtools/gffread/RSeQC/Qualimap/MultiQC 步骤 |
| `@count-...` | `rnaseq-count-flow` 内部 samtools/featureCounts/MultiQC 步骤 |
| `@denovo-assembly-...` | `rnaseq-denovo-assembly-flow` 内部 FastQC/fastp/Trinity/rnaSPAdes/seqkit/BUSCO/MultiQC 步骤 |
| `@denovo-expression-...` | `rnaseq-denovo-expression-flow` 内部 seqkit/QC/修剪/Salmon/MultiQC 步骤 |
| `@denovo-annotation-...` | `rnaseq-denovo-annotation-flow` 内部 seqkit/TransDecoder/DIAMOND 步骤 |
| `@de-...` | `rnaseq-de-flow` 内部 DESeq2/PCA/绘图步骤 |
| `@enrichment-...` | `rnaseq-enrichment-flow` 内部 ORA/GSEA/绘图步骤 |

例如，大型无参组装如果需要限制 Trinity 的 in silico read normalization 覆盖深度，
以降低临时磁盘和内存压力，可以写：

```sh
taf-rnaseq-standard-flow \
  --mode denovo \
  --samples samples.tsv \
  --metadata metadata.tsv \
  --design '~ condition' \
  --contrast condition:treated:control \
  --protein-db proteins.faa \
  --go-map protein_go_map.tsv \
  --outdir rnaseq-denovo-standard-out \
  @denovo-assembly-trinity-assembly-step: --normalize_max_read_cov 50 @:
```

这类参数适合高级用户处理资源、算法或特殊工具选项；常规分析仍应优先使用顶层参数。

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
    denovo_assembly/
    denovo_expression/
    denovo_annotation/
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
- 如果没有可靠参考基因组/注释，明确使用 `--mode denovo`，并认真准备本地 `--protein-db` 与 `--go-map`。
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

### 为什么无参模式不自动开启

缺少 `--genome` 可能是用户忘记传参，也可能是真正没有参考基因组。自动切换会把错误输入悄悄变成另一种分析路线，所以 `0.3.0-r1` 要求显式写 `--mode denovo`。

### 无参结果是不是 gene-level 结果

默认不是。无参路线的主要 feature 是 assembled transcript ID。只有当后续提供可靠 transcript-to-gene 或 clustering 映射时，才可以构建更接近 gene/pseudo-gene 层面的矩阵。`0.3.0-r1` 的 standard denovo 路线以 transcript-level DE 和 homology-derived enrichment 为主。

## 11. 边界

`rnaseq-standard-flow` 不做这些事情：

- 不联网下载 reference、annotation、gene sets、protein database、GO map 或其他数据库。
- 不自动选择生物学 design。
- 不自动判断 strandedness。
- 不自动判断项目应该走 reference 还是 denovo。
- 不自动删除失败样本。
- 不替代完整生产级 workflow engine。
- `0.3.0-r1` 可以构建 Kallisto index，但 reference 表达定量仍走 Salmon。
- 无参同源注释是功能证据，不等同于人工审定注释；assembled transcript ID 也不等同于已知 reference gene ID。

如果项目需要复杂 batch correction、多因素交互设计、splicing analysis、fusion detection、allele-specific expression、single-cell RNA-seq 或临床级报告，应在本流程输出的基础上设计额外分析。
