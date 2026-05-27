# rnaseq-standard-flow User Manual

`rnaseq-standard-flow` is the TAFFISH RNA-seq standard umbrella flow. It composes published, version-pinned RNA-seq subflows into one reproducible command that starts from common bulk RNA-seq inputs and produces expression matrices, differential expression results, enrichment results, plots, logs, version records, and a static HTML project report.

This manual is written for:

- New users who want to run a standard RNA-seq route.
- Bioinformatics users who need to prepare the input files correctly.
- Project maintainers who need auditable commands, versions, and deliverables.

## 1. What This Flow Does

The default route is Salmon-first:

```text
FASTQ + genome + annotation + metadata + GMT
-> rnaseq-index-flow
-> rnaseq-expression-flow
-> rnaseq-de-flow
-> rnaseq-enrichment-flow
-> rnaseq-report-flow
```

With `--route both`, the flow also runs an alignment/count/QC branch:

```text
rnaseq-index-flow --genome-indexer hisat2
-> rnaseq-alignment-flow
-> rnaseq-alignment-qc-flow
-> rnaseq-count-flow
```

By default, differential expression uses the Salmon/tximport gene-level count matrix. With `--route both --de-source featurecounts`, differential expression uses the featureCounts gene count matrix from the alignment/count branch.

This flow is not a Nextflow- or Snakemake-style workflow engine. It is a TAFFISH shell-native umbrella flow: each core step is executed through an explicit, version-pinned `taf-rnaseq-*-flow` dependency, and commands, versions, logs, and results are written under one explicit output directory.

## 2. Quick Start

Most projects start with:

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

To add HISAT2 BAM files, BAM-level QC, and featureCounts counts while keeping Salmon as the DE source:

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

To run DESeq2 from featureCounts counts:

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

## 3. Preparing Inputs

### 3.1 FASTQ Sample Table

`--samples` points to a tab-separated file. Minimal single-end example:

```text
sample_id	read1	condition
S1	data/S1.fastq.gz	control
S2	data/S2.fastq.gz	treated
```

Minimal paired-end example:

```text
sample_id	read1	read2	condition
S1	data/S1_R1.fastq.gz	data/S1_R2.fastq.gz	control
S2	data/S2_R1.fastq.gz	data/S2_R2.fastq.gz	treated
```

Rules:

- `sample_id` must be unique and non-empty.
- `read1` is required.
- `read2` enables paired-end mode.
- Relative FASTQ paths are resolved relative to the directory containing `samples.tsv`.
- The default Salmon route can handle either single-end or paired-end data.
- `--route both` requires a homogeneous layout: all samples must be single-end or all samples must be paired-end.

Additional columns such as `condition`, `batch`, `library_layout`, or `strandedness` are useful for review, but the formal DE design comes from `--metadata`.

### 3.2 Reference Genome FASTA

`--genome` points to a reference genome FASTA:

```text
>chrI
ACGT...
>chrII
ACGT...
```

Important points:

- Annotation sequence IDs must match the first token of each FASTA header.
- If a FASTA header is `>chrI some description`, the annotation should use `chrI`.
- For formal projects, use genome and annotation files from the same database release.

### 3.3 GFF3/GTF Annotation

`--annotation` accepts GFF3 or GTF. `rnaseq-index-flow` standardizes it and uses it to:

- extract transcript FASTA;
- write `tx2gene.tsv`;
- write gene-level GTF/GFF3 files;
- provide annotation for featureCounts and RNA-seq QC when `--route both` is used.

A common failure is inconsistent sequence naming between FASTA and annotation files, for example `chrI` in FASTA but `I` in GFF3.

### 3.4 Metadata Table

`--metadata` is passed to `rnaseq-de-flow`. The default sample column is `sample`:

```text
sample	condition	batch
S1	control	B1
S2	control	B1
S3	treated	B2
S4	treated	B2
```

Rules:

- The default `sample` column must match the sample names in the count matrix.
- Use `--sample-column` if your sample column has a different name.
- Columns such as `condition` and `batch` can be used in the DESeq2 design.
- Do not include samples that are absent from the count matrix, and do not omit count-matrix samples.

### 3.5 Design and Contrast

`--design` is a DESeq2 design formula:

```text
--design '~ condition'
```

With batch:

```text
--design '~ batch + condition'
```

`--contrast` uses `FACTOR:NUMERATOR:DENOMINATOR`:

```text
--contrast condition:treated:control
```

This means:

```text
treated / control
```

Positive log2 fold change means higher expression in treated samples. Negative log2 fold change means higher expression in control samples.

### 3.6 GMT Gene Sets

`--gene-sets` points to an offline GMT file:

```text
term_id	description	gene1	gene2	gene3
```

The flow does not download GO, KEGG, Reactome, or any other database at runtime. Prepare the GMT file ahead of time, and make sure its gene IDs use the same ID space as the DE results.

### 3.7 Background Gene Universe

`--background` is optional but recommended. A common format is:

```text
gene_id
YAL001C
YAL002W
```

When a background is provided, standard-flow filters the DE significant-gene list and ranked-gene table to that ID space before calling enrichment. It records the filtering summary in:

```text
04_reports/enrichment_background_filter.tsv
```

This helps keep enrichment tests consistent when annotation, DE results, and gene sets do not contain exactly the same set of gene IDs.

## 4. Common Run Modes

### 4.1 Default Salmon-First Route

Use this route for:

- common bulk RNA-seq quantification;
- teaching and lightweight projects;
- projects that do not need BAM deliverables;
- faster runs with a smaller resource footprint.

The default route outputs gene-level counts, TPMs, DESeq2 results, ORA/GSEA results, and the final HTML report.

### 4.2 `--route both`

Use this route when you need:

- sorted BAM files;
- HISAT2 alignment summaries;
- RSeQC, Qualimap, and SAMtools BAM-level QC;
- a traditional alignment-count evidence branch.

This route uses more CPU time, disk, and memory. DE still uses Salmon counts unless `--de-source featurecounts` is also set.

### 4.3 `--route both --de-source featurecounts`

Use this route when:

- a project requires classical genome alignment + featureCounts + DESeq2;
- differential expression must be based on BAM/featureCounts;
- you want Salmon and featureCounts as two auditable evidence lines.

This does not mean featureCounts is universally more correct than Salmon. Salmon-first and alignment-count are different modeling routes. Choose based on project goals, sample quality, annotation quality, and delivery requirements.

## 5. Parameter Reference

### 5.1 Required Parameters

| Parameter | Meaning |
| --- | --- |
| `--samples` | FASTQ sample table. |
| `--genome` | Reference genome FASTA. |
| `--annotation` | GFF3/GTF annotation. |
| `--metadata` | DESeq2 metadata table. |
| `--design` | DESeq2 design formula. |
| `--contrast` | Contrast in `factor:numerator:denominator` form. |
| `--gene-sets` | Offline GMT gene-set file. |
| `--outdir` | Dedicated output directory. Existing directories are refused by default. |

### 5.2 Main Flow Parameters

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--threads` | `2` | Threads for index, expression, and related subflows. |
| `--route` | `salmon` | `salmon` or `both`. |
| `--de-source` | `salmon` | `salmon` or `featurecounts`. |
| `--project-name` | `RNA-seq project` | Project name shown in the final HTML report. |
| `--force` | off | Replace existing standard-flow outputs under `--outdir`. |

### 5.3 Expression Parameters

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--library-type` | `A` | Salmon library type. `A` lets Salmon infer the type. |
| `--indexer` | `salmon` | Build `salmon` or `both` reference indexes. r1 quantification still uses Salmon. |
| `--kmer` | `31` | k-mer setting for Salmon/Kallisto index construction. |
| `--trim` | off | Run fastp trimming inside the expression subflow. |
| `--skip-fastqc` | off | Skip FastQC inside the expression subflow. |
| `--min-assigned-frags` | `10` | Minimum assigned fragments check for Salmon quantification. |
| `--counts-from-abundance` | `no` | tximport gene-count handling: `no`, `scaledTPM`, `lengthScaledTPM`, or `dtuScaledTPM`. |

### 5.4 Alignment, Count, and QC Parameters

These parameters are mainly used with `--route both`.

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--rna-strandness` | `none` | HISAT2 RNA strandness: `none`, `F`, `R`, `FR`, or `RF`. |
| `--alignment-min-mapq` | `0` | MAPQ filter used by alignment/count steps. |
| `--count-strand` | `0` | featureCounts strand mode: `0` unstranded, `1` stranded, `2` reversely stranded. |
| `--count-feature-type` | `exon` | Feature type counted by featureCounts. |
| `--count-attribute` | `gene_id` | Attribute used to group features. |
| `--count-min-assigned-reads` | `0` | Minimum assigned reads check for count-flow. |
| `--qc-mapq` | `30` | MAPQ cutoff used by RSeQC. |
| `--infer-sample-size` | `200000` | Sampling size for RSeQC infer_experiment.py. |
| `--java-mem-size` | `4G` | Qualimap Java memory setting. |
| `--sequencing-protocol` | `non-strand-specific` | Qualimap protocol: `non-strand-specific`, `strand-specific-forward`, or `strand-specific-reverse`. |

### 5.5 Differential Expression Parameters

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--sample-column` | `sample` | Sample column in metadata. |
| `--gene-column` | `gene_id` | Gene column in the count matrix. |
| `--padj-cutoff` | `0.05` | Adjusted P-value cutoff for significant genes. |
| `--lfc-cutoff` | `1` | Absolute log2 fold-change cutoff for significant genes. |
| `--fit-type` | `parametric` | DESeq2 dispersion fit: `parametric`, `local`, or `mean`. |
| `--lfc-shrink` | `none` | LFC shrinkage: `none`, `ashr`, or `apeglm`. |
| `--coef` | empty | Required when `--lfc-shrink apeglm` is used. |
| `--min-count` | `1` | Minimum count for low-expression filtering. |
| `--min-samples` | `2` | Minimum number of samples passing `--min-count`. |
| `--top-var` | `500` | Number of top variable genes used for sample-structure plots. |
| `--top-heatmap` | `50` | Number of top genes shown in the heatmap. |

### 5.6 Enrichment Parameters

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--enrichment-min-size` | `2` | Minimum GMT gene-set size. |
| `--enrichment-max-size` | `500` | Maximum GMT gene-set size. |
| `--enrichment-pvalue-cutoff` | `1` | Enrichment P-value cutoff. The default keeps broad output for reporting and downstream filtering. |
| `--enrichment-padj-method` | `BH` | Multiple-testing adjustment method. |
| `--enrichment-top-n` | `20` | Number of top gene sets used in enrichment plots. |
| `--enrichment-seed` | `1` | Random seed used by enrichment steps such as GSEA. |

## 6. Output Layout and Key Files

Standard output tree:

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
    plot_files.tsv
    commands.sh
    versions.tsv
    methods.txt
    flow_summary.tsv
    subflows.tsv
  run.manifest.json
```

Important files:

| File | Purpose |
| --- | --- |
| `04_reports/rnaseq_report.html` | Main project report entry point. |
| `03_results/report/` | Full report-flow collected output. |
| `03_results/plots/png/` | Top-level PNG plot collection. |
| `03_results/plots/pdf/` | Top-level PDF plot collection. |
| `04_reports/plot_files.tsv` | Top-level plot index. |
| `04_reports/commands.sh` | Command record for the umbrella flow and upstream subflows. |
| `04_reports/versions.tsv` | Version record. |
| `04_reports/methods.txt` | Methods text. |
| `04_reports/subflows.tsv` | Status of each subflow. |
| `run.manifest.json` | Machine-readable run manifest. |

## 7. Reading the Final HTML Report

The final report is a static HTML file that can be opened directly in a browser. It includes:

- one-click English/Chinese switching;
- sidebar navigation with active-section highlighting;
- project overview and key metrics;
- RNA-seq workflow diagrams;
- sections for reference, expression, alignment/count/QC, DE, enrichment, deliverables, and provenance;
- PNG/PDF plot cards;
- links to MultiQC, FastQC, Qualimap, and other HTML subreports;
- tool and TAFFISH source links.

Common plot meanings:

| Plot | Meaning |
| --- | --- |
| PCA plot | Whether samples cluster by condition, batch, or possible outlier structure. |
| Sample correlation heatmap | Expression similarity between samples. |
| Expression distribution | Whether expression distributions show strong shifts. |
| MA plot | Relationship between mean expression and log2 fold change. |
| Volcano plot | Joint view of effect size and statistical significance. |
| DEG counts barplot | Number of up- and down-regulated genes. |
| Heatmap | Expression pattern of top differential genes across samples. |
| Top genes expression | Expression values for representative differential genes. |
| ORA dotplot/barplot | Over-represented gene sets among significant genes. |
| GSEA NES plot | Directional gene-set signal in the ranked gene list. |
| GSEA enrichment curves | Running enrichment curves for representative gene sets. |

The report helps organize the analysis but does not replace project-specific biological interpretation. Final conclusions should account for experimental design, sample quality, batch structure, validation experiments, and domain knowledge.

## 8. Running With the Yeast Test Data

Use `taf-rnaseq-yeast-get-data` to prepare the yeast SNF2 example dataset in a local directory of your choice. Full preparation downloads public ENA, SGD, and GO resources, so check network access and disk space first.

```sh
taf-rnaseq-yeast-get-data \
  --outdir yeast-snf2-data-v1 \
  --stage all \
  --resume true
```

After the data flow finishes, use `<outdir>/03_results` as the data root. The example below assumes this directory exists:

```text
yeast-snf2-data-v1/03_results/
```

Set a shell variable to keep the command readable:

```sh
DATA="$PWD/yeast-snf2-data-v1/03_results"
```

Then run standard-flow:

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

Use `--stage plan` first if you only want to inspect the acquisition plan. Use `--stage reference` or `--stage genesets` if you only need those packages. In real projects, `metadata.tsv` must match the selected sample set exactly; if you manually subset samples, subset metadata at the same time.

## 9. Recommended Choices

For common projects:

- Start with the default Salmon route to confirm that FASTQ, reference, metadata, and GMT inputs line up.
- Add `--route both` when BAM files or alignment QC are required.
- Use `--route both --de-source featurecounts` when a classical count route is required for delivery.
- Provide `--background` when a reliable gene universe is available.
- Keep `--library-type A` if you are not certain about Salmon library type.
- Skip `--trim` if reads have already been cleaned by a trusted upstream process.
- Do not use very small `--top-var` or `--top-heatmap` values for formal reports unless you are only running a quick test.

## 10. Troubleshooting

### Annotation and genome sequence IDs do not match

Symptom: index-flow fails during reference validation or transcript extraction.

Fix: make the first token of each FASTA header match the sequence IDs in GFF3/GTF.

### Metadata samples do not match the count matrix

Symptom: DE flow reports a mismatch between metadata and count matrix samples.

Fix: check `sample_id`, the metadata `sample` column, and `--sample-column`.

### GMT, background, and DE gene IDs are inconsistent

Symptom: enrichment results are sparse or empty, or background filtering retains very few genes.

Fix: make sure the count matrix, annotation, GMT file, and background file use the same gene ID system.

### Mixed single-end and paired-end samples

Symptom: `--route both` is refused.

Fix: split the data into separate projects or normalize the layout during data preparation.

### Output directory already exists

Symptom: the flow refuses an existing `--outdir`.

Fix: use a new output directory, or use `--force` only after confirming that replacement is intended.

### Report links do not open

Open the final report from `04_reports/rnaseq_report.html`, and keep the whole output tree together. The HTML subreports and plot links depend on the surrounding output directory structure.

## 11. Boundaries

`rnaseq-standard-flow` does not:

- download references, annotations, gene sets, or databases at runtime;
- choose the biological design automatically;
- infer strandedness automatically;
- remove failed samples automatically;
- replace a production workflow engine;
- use Kallisto for expression quantification in r1, although it can build a Kallisto index.

If a project needs complex batch correction, interaction models, splicing analysis, fusion detection, allele-specific expression, single-cell RNA-seq, or clinical reporting, use this flow as a reproducible baseline and design additional analyses on top of its outputs.
