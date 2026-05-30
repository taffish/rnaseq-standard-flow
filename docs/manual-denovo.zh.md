# rnaseq-standard-flow 无参分析手册

本手册专门说明 `rnaseq-standard-flow 0.2.0-r1` 的无参路线，也就是
`--mode denovo`。它面向没有可靠参考基因组、没有高质量基因注释，或者研究问题更偏向
转录本发现的 bulk RNA-seq 项目。

完整 standard-flow 总手册见 `docs/manual.zh.md`。

## 1. 什么时候使用无参模式

当研究对象不能被可靠的 reference genome 和 gene annotation 表达时，可以考虑
`--mode denovo`。常见场景包括非模式物种、公共参考质量很差、近缘物种差异较大、基因注释不完整，
或者项目目标本来就是发现转录本。

不要只是因为临时找不到 `--genome` 或 `--annotation` 就默认使用无参。缺少参考输入可能是用户漏传文件，
也可能是真的没有参考。为了避免把错误输入静默变成另一种分析策略，流程不会自动切换到无参，必须显式写：

```sh
--mode denovo
```

无参路线回答的问题是：

```text
哪些组装出来的 transcript 在不同条件之间表达发生变化？
这些 transcript 可以通过同源蛋白和 GO 证据关联到哪些功能过程？
```

它默认不能直接回答：

```text
哪些已经审定的 reference gene 发生变化？
```

除非后续提供可靠的 transcript-to-gene 或 transcript clustering 映射。

## 2. 生物学和技术路线

无参 standard route 是：

```text
FASTQ + metadata + protein FASTA + protein-to-GO map
-> rnaseq-denovo-assembly-flow
-> rnaseq-denovo-expression-flow
-> rnaseq-denovo-annotation-flow
-> rnaseq-de-flow
-> rnaseq-enrichment-flow
-> rnaseq-report-flow
```

每一步的生物学含义是：

- assembly 在没有 genome 的情况下，从 reads 重建 transcript sequence。
- expression quantification 把 reads 定量到 assembled transcripts 上。
- ORF prediction 和 homology annotation 把 transcript 推断出的蛋白与已知蛋白证据连接起来。
- differential expression 比较不同条件下 transcript-level counts 的变化。
- enrichment 检测发生变化的 transcript 所关联的 GO term 是否富集。
- final report 会把无参证据与 reference gene-level 证据区分开。

这条路线有价值，但解释时比有参路线更依赖证据质量。assembly 质量、测序深度、样本覆盖、
transcript fragmentation、isoform redundancy、蛋白数据库选择和 GO 映射质量都会影响最终结论。

## 3. 最小命令

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

正式数据通常还应设置线程、内存和项目名：

```sh
taf-rnaseq-standard-flow \
  --mode denovo \
  --samples samples.tsv \
  --metadata metadata.tsv \
  --design '~ batch + condition' \
  --contrast condition:treated:control \
  --protein-db proteins.faa \
  --go-map protein_go_map.tsv \
  --threads 16 \
  --max-memory 64G \
  --project-name "No-reference RNA-seq" \
  --outdir rnaseq-denovo-standard-out
```

## 4. 必需输入

### 4.1 FASTQ 样本表

`--samples samples.tsv` 是 tab 分隔表，每行一个生物样本。

必需列：

```text
sample_id	read1
```

双端测序可以增加 `read2`：

```text
sample_id	read1	read2
WT_01	reads/WT_01_R1.fq.gz	reads/WT_01_R2.fq.gz
WT_02	reads/WT_02_R1.fq.gz	reads/WT_02_R2.fq.gz
KO_01	reads/KO_01_R1.fq.gz	reads/KO_01_R2.fq.gz
KO_02	reads/KO_02_R1.fq.gz	reads/KO_02_R2.fq.gz
```

同一次运行不能混合 single-end 和 paired-end 样本。相对 FASTQ 路径会基于
`samples.tsv` 所在目录解析。

对于无参组装，reads 的平衡性和代表性很重要。如果某个条件测序深度明显更高，assembler
可能更容易恢复该条件中的 transcript。流程通过 `--no-normalize` 暴露 Trinity 的
normalization 控制，但默认保留 Trinity 的 normalization 行为。

### 4.2 Metadata 表

`--metadata metadata.tsv` 被 `rnaseq-de-flow` 使用，必须包含样本 ID 列。默认列名是
`sample`。

示例：

```text
sample	condition	batch
WT_01	WT	B1
WT_02	WT	B2
KO_01	KO	B1
KO_02	KO	B2
```

metadata 中的样本名必须和 `samples.tsv` 的 `sample_id` 一致。如果手动筛选 FASTQ
样本，也必须同步筛选 metadata。

### 4.3 Design 和 Contrast

`--design` 是 DESeq2 design formula。最简单写法：

```sh
--design '~ condition'
```

如果实验有已知 batch，可以写：

```sh
--design '~ batch + condition'
```

`--contrast` 指定比较：

```sh
--contrast condition:treated:control
```

其中 numerator 相对 denominator 的 log2 fold-change 为正。

### 4.4 蛋白数据库

`--protein-db proteins.faa` 是本地蛋白 FASTA，供 `rnaseq-denovo-annotation-flow`
通过 DIAMOND 做同源证据搜索。

数据库选择建议：

- 非模式物种优先使用近缘、注释质量高的物种或类群蛋白集。
- 酵母类数据优先使用高质量 yeast protein set，而不是无关的大数据库。
- 探索性项目可以使用更大的数据库提高敏感性，但运行时间和歧义也会增加。

FASTA header ID 必须和 `--go-map` 中的 `subject_id` 对得上。例如：

```text
>P12345 hypothetical kinase
MSTNPKPQR...
```

那么 `protein_go_map.tsv` 中应使用 `P12345` 作为 `subject_id`。

### 4.5 Protein-to-GO Map

`--go-map protein_go_map.tsv` 把 protein subject ID 映射到 GO term。推荐列：

```text
subject_id	go_id	go_name	namespace
P12345	GO:0004672	protein kinase activity	molecular_function
P12345	GO:0006468	protein phosphorylation	biological_process
P67890	GO:0005737	cytoplasm	cellular_component
```

必需列：

- `subject_id`
- `go_id`

推荐可选列：

- `go_name`
- `namespace`

denovo annotation flow 会用 DIAMOND hit 加上这个 GO map 生成：

```text
03_results/denovo_annotation/03_results/gene_sets/denovo_go.gmt
03_results/denovo_annotation/03_results/gene_sets/denovo_background.tsv
```

这些文件随后进入 enrichment。若 GO map 覆盖很低，富集结果会稀疏或偏倚。

## 5. 无参参数

| 参数 | 默认值 | 含义 |
| --- | --- | --- |
| `--mode` | `reference` | 无参分析必须设置为 `denovo`。 |
| `--assembler` | `trinity` | `trinity` 或 `rnaspades`；传给 `rnaseq-denovo-assembly-flow`。 |
| `--max-memory` | `4G` | assembler 内存限制；正式数据通常需要增大。 |
| `--min-contig-len` | `200` | 保留 assembled transcript 的最小长度。 |
| `--ss-lib-type` | `none` | strand-specific library type，传给支持该参数的 assembly 工具。 |
| `--no-normalize` | off | 添加 Trinity `--no_normalize_reads`；主要用于 tiny test 或已预处理数据。 |
| `--denovo-min-orf-aa` | `50` | TransDecoder ORF 最小氨基酸长度。 |
| `--denovo-evalue` | `1e-5` | DIAMOND blastp e-value 阈值。 |
| `--denovo-max-target-seqs` | `1` | 每个预测蛋白保留的 DIAMOND subject hit 数。 |

实际建议：

- 默认优先使用 `--assembler trinity`，除非有明确理由选择 rnaSPAdes。
- 根据样本数、reads 深度和物种复杂度提高 `--max-memory`。
- `--min-contig-len 200` 适合一般项目；提高该值可减少短片段，但可能丢失短转录本或 partial transcript。
- 只有明确知道建库方向时才设置 `--ss-lib-type`。
- `--denovo-max-target-seqs 1` 可以让 enrichment 资源更清晰；只有需要更广泛注释证据且能处理歧义时才增大。

## 6. 与有参模式共享的参数

`--threads` 控制子流程线程数。无参组装通常较耗资源，正式项目一般需要高于默认值 `2`。

`--trim` 会在无参子流程中启用 fastp trimming。如果 reads 仍有 adapter 或低质量尾部，建议启用。
如果上游已经完成可靠清洗，可以不启用。

`--skip-fastqc` 会跳过 FastQC。第一次项目分析通常不建议跳过；如果上游 QC 已经归档，可以考虑。

`--padj-cutoff`、`--lfc-cutoff`、`--fit-type`、`--lfc-shrink`、`--min-count`、
`--min-samples`、`--top-var` 和 `--top-heatmap` 会传给 `rnaseq-de-flow`。
统计含义与有参模式一致，但 feature 是 assembled transcript。

`--enrichment-min-size`、`--enrichment-max-size`、`--enrichment-pvalue-cutoff`、
`--enrichment-padj-method`、`--enrichment-top-n` 和 `--enrichment-seed`
会传给 `rnaseq-enrichment-flow`。无参模式中的 GMT/background 来自同源注释。

## 7. 无参模式不使用的参数

这些参数只用于 reference 模式：

- `--genome`
- `--annotation`
- `--gene-sets`
- `--background`
- `--route both`
- `--de-source featurecounts`
- `--rna-strandness`
- `--alignment-min-mapq`
- `--count-strand`
- `--count-feature-type`
- `--count-attribute`
- `--count-min-assigned-reads`
- `--qc-mapq`
- `--infer-sample-size`
- `--java-mem-size`
- `--sequencing-protocol`

原因既是技术上的，也是生物学上的：HISAT2 alignment、featureCounts gene counting、
RSeQC 和 Qualimap 都需要 reference genome 和 annotation。无参路线使用 assembled transcript
sequence 和 transcript-level matrix。

## 8. 使用 yeast 测试数据跑无参路线

`rnaseq-yeast-get-data` 获取的 yeast SNF2 数据可以同时用于有参和无参示例。推荐先按有参手册相同方式准备数据：

```sh
taf-rnaseq-yeast-get-data \
  --outdir yeast-snf2-data-v1 \
  --stage all \
  --resume true

DATA="$PWD/yeast-snf2-data-v1/03_results"
```

这样只需要获取一次数据：

- 有参示例使用同一 `DATA` 下的 genome、annotation、gene sets 和 FASTQ。
- 无参示例不把 genome 或 annotation 传给流程，只复用 FASTQ，并从同一 SGD source tarball 派生测试用 protein DB 和 GO map。

注意：yeast 有高质量参考基因组，真实 yeast 项目通常优先使用有参路线。这里的无参示例主要用于验证无参流程的工程接口、
报告结构和结果解释方式。

为了让本地测试更快，可以抽取 2 个 WT 和 2 个 SNF2KO 样本，并截取每个样本前 50,000 条 reads：

```sh
WORK="$PWD/yeast-denovo-test-inputs"
mkdir -p "$WORK/reads" "$WORK/resources"

{
  printf 'sample_id\tread1\tcondition\n'
  awk -F '\t' -v OFS='\t' '
    NR == 1 {
      for (i = 1; i <= NF; i++) col[$i] = i
      next
    }
    $(col["condition"]) == "snf2_KO" && ko < 2 {
      print $(col["sample_id"]), $(col["read1"]), $(col["condition"])
      ko++
      next
    }
    $(col["condition"]) == "WT" && wt < 2 {
      print $(col["sample_id"]), $(col["read1"]), $(col["condition"])
      wt++
      next
    }
  ' "$DATA/yeast-snf2-fastq-mini-v1/samples.tsv"
} > "$WORK/selected.tsv"

printf 'sample_id\tread1\tcondition\n' > "$WORK/samples.tsv"
printf 'sample\tcondition\n' > "$WORK/metadata.tsv"

tail -n +2 "$WORK/selected.tsv" |
while IFS="$(printf '\t')" read -r sample_id read1 condition; do
  case "$read1" in
    /*) fq="$read1" ;;
    *) fq="$DATA/yeast-snf2-fastq-mini-v1/$read1" ;;
  esac
  gzip -cd "$fq" \
    | awk 'NR <= 200000 { print }' > "$WORK/reads/$sample_id.fq"
  printf '%s\t%s\t%s\n' "$sample_id" "reads/$sample_id.fq" "$condition" >> "$WORK/samples.tsv"
  printf '%s\t%s\n' "$sample_id" "$condition" >> "$WORK/metadata.tsv"
done
```

然后从同一数据包中的 SGD reference source 派生无参注释所需的 protein FASTA 和 protein-to-GO map：

```sh
REFERENCE_TAR="$DATA/yeast-reference-sgd-r64.4.1-v1/source/S288C_reference_genome_R64-4-1_20230830.tgz"
GO_TERMS="$DATA/yeast-sgd-go-gene-sets-r64.4.1-v1/metadata/go_terms.tsv"
PROTEIN_DB="$WORK/resources/yeast_orf_trans_all_R64-4-1.protein.faa"
GO_MAP="$WORK/resources/yeast_sgd_go_map.tsv"

tar -xOzf "$REFERENCE_TAR" \
  S288C_reference_genome_R64-4-1_20230830/orf_trans_all_R64-4-1_20230830.fasta.gz \
  | gzip -cd \
  | awk '
      /^>/ { print; next }
      {
        gsub(/\*/, "")
        gsub(/[[:space:]]/, "")
        if ($0 != "") print
      }
    ' > "$PROTEIN_DB"

tar -xOzf "$REFERENCE_TAR" \
  S288C_reference_genome_R64-4-1_20230830/gene_association_R64-4-1_20230830.sgd.gz \
  | gzip -cd \
  | awk -F '\t' -v OFS='\t' -v terms="$GO_TERMS" '
      BEGIN {
        while ((getline line < terms) > 0) {
          n = split(line, t, "\t")
          if (n < 3 || t[1] == "go_id") continue
          go_name[t[1]] = t[2]
          go_namespace[t[1]] = t[3]
        }
        close(terms)
        print "subject_id", "go_id", "go_name", "namespace"
      }
      /^!/ || NF < 11 { next }
      {
        split($11, synonyms, "|")
        subject = synonyms[1]
        go = $5
        if (subject == "" || go == "") next
        name = (go in go_name) ? go_name[go] : go
        namespace = (go in go_namespace) ? go_namespace[go] : $9
        if (namespace == "P") namespace = "biological_process"
        else if (namespace == "F") namespace = "molecular_function"
        else if (namespace == "C") namespace = "cellular_component"
        key = subject SUBSEP go
        if (!(key in seen)) {
          seen[key] = 1
          print subject, go, name, namespace
        }
      }
    ' > "$GO_MAP"
```

运行无参 standard-flow：

```sh
taf-rnaseq-standard-flow \
  --mode denovo \
  --samples "$WORK/samples.tsv" \
  --metadata "$WORK/metadata.tsv" \
  --design '~ condition' \
  --contrast condition:snf2_KO:WT \
  --protein-db "$PROTEIN_DB" \
  --go-map "$GO_MAP" \
  --outdir yeast-denovo-standard-out \
  --threads 4 \
  --max-memory 8G \
  --min-contig-len 100 \
  --denovo-min-orf-aa 30 \
  --denovo-evalue 1e-3 \
  --denovo-max-target-seqs 1 \
  --no-normalize \
  --project-name "Yeast SNF2 RNA-seq de novo standard test"
```

完成后重点查看：

```text
yeast-denovo-standard-out/03_results/denovo_assembly/
yeast-denovo-standard-out/03_results/denovo_expression/
yeast-denovo-standard-out/03_results/denovo_annotation/
yeast-denovo-standard-out/03_results/de/
yeast-denovo-standard-out/03_results/enrichment/
yeast-denovo-standard-out/04_reports/rnaseq_report.html
yeast-denovo-standard-out/04_reports/report_interpretation.html
```

如果要用全量 24 个样本和每样本 500k reads 测试，可以直接基于
`$DATA/yeast-snf2-fastq-mini-v1/samples.tsv` 和
`$DATA/yeast-snf2-fastq-mini-v1/metadata.tsv` 运行；但无参组装会明显更耗时、耗内存。

## 9. 输出目录

重要目录：

```text
rnaseq-denovo-standard-out/
  00_inputs/
  01_logs/
  03_results/
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
    flow_summary.tsv
    subflows.tsv
    plot_files.tsv
    commands.sh
    versions.tsv
    methods.txt
  run.manifest.json
```

关键生物学结果：

- `denovo_assembly` 中的 assembled transcriptome FASTA。
- `denovo_expression` 中的 transcript count 和 TPM 矩阵。
- `denovo_annotation` 中的 ORF、DIAMOND 和 GO-derived annotation。
- `de` 中的 transcript-level DESeq2 表格和图。
- `enrichment` 中的 ORA/GSEA 富集表和图。
- `04_reports/rnaseq_report.html` 中的双语静态项目报告。

## 10. 如何解读结果

count matrix 是 transcript-level。每个 row ID 是 assembled transcript feature，
不一定是人工审定的 gene。多个 transcript 可能对应 isoform、paralog、fragmented piece、
allelic variant、assembly artifact 或 redundant assembly。

DE 表应解读为：

```text
assembled transcript X 在两个条件之间 read support 上升或下降
```

不要直接解读为：

```text
已知 gene Y 发生变化
```

annotation 表把 transcript 连接到 predicted ORF 和 protein hit。这是证据转移，不是人工审定。
更可靠的解释通常来自较长 ORF、质量较好的 DIAMOND hit、一致的 GO term，以及生物学重复之间稳定的信号。

enrichment 依赖生成的 denovo background。如果许多 transcript 没有 protein 或 GO 支持，
被检验的 universe 实际上只是可注释的那一部分 assembly。

## 11. 质量检查

运行后建议检查：

- assembly summary：transcript 总数、长度分布、N50 类指标。
- assembly flow 中可用的 BUSCO 输出。
- Salmon quantification 日志和 assignment 情况。
- DESeq2 的 PCA 和 sample correlation 图。
- significant transcripts 数量，以及是否被单个样本主导。
- enrichment background size 和 GO term 覆盖情况。
- `04_reports/versions.tsv` 和 `04_reports/commands.sh` 的溯源记录。

如果短 transcript 很多、ORF recovery 很弱、protein hit 覆盖差，或 GO map 极度稀疏，
下游生物学解释应更保守。

## 12. 常见问题

### 流程拒绝 `--route both`

`--route both` 只适用于 reference 模式。如果需要 genome alignment 和 featureCounts 证据，
应使用有参分析。

### 提示缺少 `--protein-db` 或 `--go-map`

无参 enrichment 来自同源证据。`--mode denovo` 必须提供这两个文件；流程不会联网下载。

### 富集结果很少

检查 DIAMOND hit 是否覆盖足够多 transcript，以及 `--go-map` 的 `subject_id`
是否能匹配 `--protein-db` 的 FASTA header ID。

### DE 结果里的 transcript ID 很难解释

这是无参分析的正常情况。应使用 annotation 输出把 transcript ID 连接到 ORF、protein hit
和 GO term。如果项目需要 gene-like ID，需要额外构建 transcript clustering 或 transcript-to-gene map。

### 组装很慢或内存不足

增大 `--max-memory`，提高 `--threads`，检查磁盘空间，并考虑先做 pilot 或 reads 归一化。
不要把 tiny test 参数用于正式数据。

## 13. 边界

无参路线不做这些事情：

- 不下载 protein database、GO map、BUSCO lineage 或其他资源。
- 不生成人工审定的 gene annotation。
- 不保证一个 assembled transcript 对应一个 biological gene。
- 不运行 reference-genome alignment 或 featureCounts。
- 不把 homology-derived annotation 等同于人工注释。
- 不自动判断项目是否应该使用无参分析。

无参模式适合作为可复现的 no-reference baseline。若用于发表级项目，应在形成最终结论前认真审查
assembly 质量、annotation 覆盖度和生物学合理性。
