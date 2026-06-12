# rnaseq-standard-flow User Manual

`rnaseq-standard-flow` is the TAFFISH RNA-seq standard umbrella flow. It composes published, version-pinned RNA-seq subflows into one reproducible command that starts from common bulk RNA-seq inputs and produces expression matrices, differential expression results, enrichment results, plots, logs, version records, and a static HTML project report.

This manual is written for:

- New users who want to run a standard RNA-seq route.
- Bioinformatics users who need to prepare the input files correctly.
- Project maintainers who need auditable commands, versions, and deliverables.

## 1. What This Flow Does

Version `0.3.0-r1` keeps the reference / Salmon-first route as the default
and remains compatible with the previous reference command line:

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

With explicit `--mode denovo`, the flow runs the no-reference route:

```text
FASTQ + metadata + local protein FASTA + protein-to-GO map
-> rnaseq-denovo-assembly-flow
-> rnaseq-denovo-expression-flow
-> rnaseq-denovo-annotation-flow
-> rnaseq-de-flow
-> rnaseq-enrichment-flow
-> rnaseq-report-flow
```

The de novo route is never triggered automatically by missing `--genome` or
`--annotation`. You must select it with `--mode denovo` to avoid silently
turning a missing reference input into a different analysis strategy.

For projects that primarily use the no-reference route, also read the
dedicated de novo manual:

- `docs/manual-denovo.en.md`

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

For a project without a reliable reference genome and annotation, use explicit
de novo mode with a local protein FASTA and protein-to-GO map:

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

### 3.2 Reference Genome FASTA (Reference Mode)

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

### 3.3 GFF3/GTF Annotation (Reference Mode)

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

### 3.6 GMT Gene Sets (Reference Mode)

`--gene-sets` points to an offline GMT file:

```text
term_id	description	gene1	gene2	gene3
```

The flow does not download GO, KEGG, Reactome, or any other database at runtime. Prepare the GMT file ahead of time, and make sure its gene IDs use the same ID space as the DE results.

### 3.7 Background Gene Universe (Reference Mode)

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

### 3.8 Protein Database and GO Map for De Novo Mode

`--mode denovo` does not use `--genome`, `--annotation`, `--gene-sets`, or
`--background`. It requires:

```text
--protein-db proteins.faa
--go-map protein_go_map.tsv
```

`--protein-db` is a local protein FASTA. `rnaseq-denovo-annotation-flow` builds
a local DIAMOND database from it and compares proteins predicted from assembled
transcripts against that database.

`--go-map` maps protein subject IDs to GO terms. Recommended format:

```text
subject_id	go_id	go_name	namespace
PROT1	GO:0006412	translation	biological_process
PROT2	GO:0005737	cytoplasm	cellular_component
```

The de novo route uses DIAMOND subject hits plus this GO map to generate:

```text
03_results/denovo_annotation/03_results/gene_sets/denovo_go.gmt
03_results/denovo_annotation/03_results/gene_sets/denovo_background.tsv
```

Those files are then passed to the enrichment subflow. The IDs are assembled
transcript IDs, not known reference gene IDs. Homology-transferred annotation
is evidence, not manually curated function.

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

### 4.4 `--mode denovo`

Use this route when:

- the organism has no reliable reference genome;
- a reference exists but gene annotation is poor and transcript discovery is important;
- you need a transcriptome-level expression matrix for a non-model organism.

The de novo route assembles transcripts first, then quantifies reads against
assembled transcripts with Salmon. Differential expression uses the
transcript-level count matrix and fixes `--gene-column transcript_id`. The
features are assembled transcripts rather than reference genes. Downstream
enrichment depends on the biological relevance and quality of `--protein-db`
and `--go-map`.

De novo mode cannot use `--route both` or `--de-source featurecounts`, because
BAM alignment, featureCounts, and alignment QC require a reference genome and
annotation.

## 5. Parameter Reference

### 5.1 Required Parameters

| Parameter | Meaning |
| --- | --- |
| `--samples` | FASTQ sample table. |
| `--mode` | `reference` or `denovo`. Default: `reference`, preserving previous reference-route compatibility. |
| `--genome` | Required in reference mode; reference genome FASTA. |
| `--annotation` | Required in reference mode; GFF3/GTF annotation. |
| `--metadata` | DESeq2 metadata table. |
| `--design` | DESeq2 design formula. |
| `--contrast` | Contrast in `factor:numerator:denominator` form. |
| `--gene-sets` | Required in reference mode; offline GMT gene-set file. |
| `--protein-db` | Required in denovo mode; local protein FASTA. |
| `--go-map` | Required in denovo mode; protein subject ID to GO term TSV. |
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
| `--indexer` | `salmon` | Build `salmon` or `both` reference indexes. Reference expression quantification still uses Salmon. |
| `--kmer` | `31` | k-mer setting for Salmon/Kallisto index construction. |
| `--trim` | off | Run fastp trimming inside the expression subflow. |
| `--skip-fastqc` | off | Skip FastQC inside the expression subflow. |
| `--min-assigned-frags` | `10` | Minimum assigned fragments check for Salmon quantification. |
| `--counts-from-abundance` | `no` | tximport gene-count handling: `no`, `scaledTPM`, `lengthScaledTPM`, or `dtuScaledTPM`. |

### 5.4 De Novo Parameters

These parameters are used only with `--mode denovo`.

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--assembler` | `trinity` | Assembler: `trinity` or `rnaspades`. |
| `--max-memory` | `4G` | Assembler memory limit, for example `4G` or `32G`. |
| `--min-contig-len` | `200` | Minimum assembled transcript length retained. |
| `--ss-lib-type` | `none` | Strand-specific library type: `none`, `F`, `R`, `FR`, or `RF`. |
| `--no-normalize` | off | Pass Trinity `--no_normalize_reads` through assembly-flow. |
| `--denovo-min-orf-aa` | `50` | Minimum TransDecoder ORF amino-acid length. |
| `--denovo-evalue` | `1e-5` | DIAMOND blastp e-value cutoff. |
| `--denovo-max-target-seqs` | `1` | Number of DIAMOND subject hits retained per predicted protein. |

For real projects, set `--max-memory` according to sequencing depth and
transcriptome complexity. `--no-normalize` is mainly for tiny tests or already
downsampled data; it is not a universal default for production projects.

### 5.5 Alignment, Count, and QC Parameters

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

### 5.6 Differential Expression Parameters

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--sample-column` | `sample` | Sample column in metadata. |
| `--gene-column` | `gene_id` | Reference-mode gene column in the count matrix; denovo mode fixes this to `transcript_id`. |
| `--padj-cutoff` | `0.05` | Adjusted P-value cutoff for significant genes. |
| `--lfc-cutoff` | `1` | Absolute log2 fold-change cutoff for significant genes. |
| `--fit-type` | `parametric` | DESeq2 dispersion fit: `parametric`, `local`, or `mean`. |
| `--lfc-shrink` | `none` | LFC shrinkage: `none`, `ashr`, or `apeglm`. |
| `--coef` | empty | Required when `--lfc-shrink apeglm` is used. |
| `--min-count` | `1` | Minimum count for low-expression filtering. |
| `--min-samples` | `2` | Minimum number of samples passing `--min-count`. |
| `--top-var` | `500` | Number of top variable genes used for sample-structure plots. |
| `--top-heatmap` | `50` | Number of top genes shown in the heatmap. |

### 5.7 Enrichment Parameters

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--enrichment-min-size` | `2` | Minimum GMT gene-set size. |
| `--enrichment-max-size` | `500` | Maximum GMT gene-set size. |
| `--enrichment-pvalue-cutoff` | `1` | Enrichment P-value cutoff. The default keeps broad output for reporting and downstream filtering. |
| `--enrichment-padj-method` | `BH` | Multiple-testing adjustment method. |
| `--enrichment-top-n` | `20` | Number of top gene sets used in enrichment plots. |
| `--enrichment-seed` | `1` | Random seed used by enrichment steps such as GSEA. |

### 5.8 Advanced Internal Tool Passthrough

Most projects do not need this section. `rnaseq-standard-flow 0.3.0-r1`
exposes routine options as stable top-level parameters. For unusual projects,
a lower-level tool inside a subflow may still need an option that is not
promoted to the top-level interface. Standard-flow now bridges those internal
subflow `@step:` blocks with namespaced top-level blocks.

The naming pattern is:

```text
@<standard-step>-<child-tool-step>: ... @:
```

For example, `@denovo-assembly-trinity-assembly-step:` is passed to
`rnaseq-denovo-assembly-flow` as its internal `@trinity-assembly-step:` block.
All blocks are empty by default and do not alter existing command behavior
unless explicitly supplied. The complete block list is in `docs/help.md` and
the README. Common groups are:

| Group prefix | Subflow controlled |
| --- | --- |
| `@index-...` | `rnaseq-index-flow` internal AGAT/gffread/indexer steps |
| `@expression-...` | `rnaseq-expression-flow` internal QC/trimming/Salmon/tximport/MultiQC steps |
| `@alignment-...` | `rnaseq-alignment-flow` internal fastp/HISAT2/samtools/MultiQC steps |
| `@alignment-qc-...` | `rnaseq-alignment-qc-flow` internal samtools/gffread/RSeQC/Qualimap/MultiQC steps |
| `@count-...` | `rnaseq-count-flow` internal samtools/featureCounts/MultiQC steps |
| `@denovo-assembly-...` | `rnaseq-denovo-assembly-flow` internal FastQC/fastp/Trinity/rnaSPAdes/seqkit/BUSCO/MultiQC steps |
| `@denovo-expression-...` | `rnaseq-denovo-expression-flow` internal seqkit/QC/trimming/Salmon/MultiQC steps |
| `@denovo-annotation-...` | `rnaseq-denovo-annotation-flow` internal seqkit/TransDecoder/DIAMOND steps |
| `@de-...` | `rnaseq-de-flow` internal DESeq2/PCA/plot steps |
| `@enrichment-...` | `rnaseq-enrichment-flow` internal ORA/GSEA/rendering steps |

For example, a large de novo assembly run may need stricter Trinity in silico
read normalization to reduce temporary disk and memory pressure:

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

Use passthrough for advanced resource, algorithm, or tool-specific tuning.
Prefer documented top-level options for routine analyses.

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

Important files:

| File | Purpose |
| --- | --- |
| `04_reports/rnaseq_report.html` | Main project report entry point. |
| `04_reports/report_interpretation.html` | Companion RNA-seq interpretation guide with long-form biological and technical explanations. |
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
- Use explicit `--mode denovo` when no reliable reference genome/annotation is available, and prepare a local `--protein-db` plus `--go-map` carefully.
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

### Why de novo mode is not automatic

Missing `--genome` may mean the user forgot an input, or it may mean the
project truly lacks a reference. Automatic switching would silently change the
analysis strategy, so `0.3.0-r1` requires explicit `--mode denovo`.

### Are de novo results gene-level results?

Not by default. The primary features are assembled transcript IDs. Gene-like
or pseudo-gene matrices require a reliable transcript-to-gene or clustering
map. The `0.3.0-r1` standard de novo route is transcript-level DE plus
homology-derived enrichment.

## 11. Boundaries

`rnaseq-standard-flow` does not:

- download references, annotations, gene sets, protein databases, GO maps, or other databases at runtime;
- choose the biological design automatically;
- infer strandedness automatically;
- infer whether a project should use reference or de novo mode;
- remove failed samples automatically;
- replace a production workflow engine;
- use Kallisto for reference-mode expression quantification in `0.3.0-r1`, although it can build a Kallisto index;
- turn assembled transcript IDs into known reference gene IDs without an explicit mapping.

De novo homology annotation is evidence, not manual curation, and assembled
transcript IDs are not the same as known reference gene IDs.

If a project needs complex batch correction, interaction models, splicing analysis, fusion detection, allele-specific expression, single-cell RNA-seq, or clinical reporting, use this flow as a reproducible baseline and design additional analyses on top of its outputs.
