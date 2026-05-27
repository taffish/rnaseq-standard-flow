# rnaseq-standard-flow

`rnaseq-standard-flow` is the r1 TAFFISH RNA-seq umbrella flow. It composes
the published RNA-seq subflows into a default Salmon-first route:

1. `rnaseq-index-flow`
2. `rnaseq-expression-flow`
3. `rnaseq-de-flow`
4. `rnaseq-enrichment-flow`
5. `rnaseq-report-flow`

It can also run the explicit alignment/count branch with `--route both`:

1. `rnaseq-index-flow --genome-indexer hisat2`
2. `rnaseq-alignment-flow`
3. `rnaseq-alignment-qc-flow`
4. `rnaseq-count-flow`

By default, DESeq2 uses the Salmon/tximport gene count matrix. With
`--route both --de-source featurecounts`, DESeq2 instead uses the featureCounts
matrix from the alignment/count branch.

The flow starts from explicit local FASTQ, genome, annotation, metadata, and
GMT gene-set inputs. It does not download reference data, gene sets, models, or
databases during normal execution.

## Manuals

- [中文用户手册](docs/manual.zh.md)
- [English user manual](docs/manual.en.md)

## Usage

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

To run the alignment/count branch while keeping Salmon as the DE source:

```sh
taf-rnaseq-standard-flow \
  --samples samples.tsv \
  --genome genome.fa \
  --annotation genes.gff3 \
  --metadata metadata.tsv \
  --design '~ condition' \
  --contrast condition:treated:control \
  --gene-sets gene_sets.gmt \
  --outdir rnaseq-standard-out \
  --route both
```

To make DESeq2 use featureCounts counts:

```sh
taf-rnaseq-standard-flow \
  --samples samples.tsv \
  --genome genome.fa \
  --annotation genes.gff3 \
  --metadata metadata.tsv \
  --design '~ condition' \
  --contrast condition:treated:control \
  --gene-sets gene_sets.gmt \
  --outdir rnaseq-standard-out \
  --route both \
  --de-source featurecounts
```

`samples.tsv` must contain `sample_id` and `read1`; optional `read2` enables
paired-end mode. Relative FASTQ paths are resolved relative to the original
`samples.tsv`, not relative to the output directory.

When `--route both` is used, all rows must be homogeneous single-end or
homogeneous paired-end. Mixed single/paired sample tables are refused because
the current alignment QC and featureCounts branch uses a global paired-mode
setting.

`metadata.tsv` is passed to `rnaseq-de-flow`. By default, its sample ID column
is `sample`; use `--sample-column` if your metadata uses another column name.
The metadata samples must match the expression matrix samples generated from
`samples.tsv`.

When `--background` is provided, the standard flow filters the DE significant
gene list and ranked-gene table to that background ID space before calling
`rnaseq-enrichment-flow`. This keeps the umbrella route compatible with the
enrichment subflow's strict background contract and records the filter summary
in `04_reports/enrichment_background_filter.tsv`.

## Outputs

The standard flow writes only under the explicit `--outdir`:

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
    enrichment_background_filter.tsv
  run.manifest.json
```

The top-level `04_reports/rnaseq_report.html` is copied from
`rnaseq-report-flow` r3; the full collector output remains available under
`03_results/report/`. The report is a branded static HTML project report with
one-click English/Chinese switching, workflow sections, plot cards, table
previews, linked QC/report HTML bundles, tool links, methods, versions, and
provenance. The r3 report adds active sidebar navigation, static workflow
diagrams, a deliverables/output-structure section, TAFFISH source links, and
separate ORA/GSEA enrichment visual summaries.

The top-level `03_results/plots/` directory is assembled by
`rnaseq-standard-flow` itself before the report collector runs. It keeps PNG
and PDF outputs in separate `png/` and `pdf/` subdirectories, while
`04_reports/plot_files.tsv` records the module, plot name, source path, copied
path, format, and file size. The top-level HTML report links back into the
self-contained `03_results/report/` collector assets, so copied report links
remain valid from `04_reports/rnaseq_report.html`.

The standard plot collection includes DE r2 plots and the enrichment r3 plot
suite: readable enrichment dotplot, classic dotplot, ORA top-term barplot,
GSEA NES plot, and GSEA enrichment curves. A full run with both PNG and PDF
formats produces 28 copied plot files.

## Dependencies

The flow depends on version-pinned TAFFISH subflows:

- `taf-rnaseq-index-flow = 0.1.0-r1`
- `taf-rnaseq-expression-flow = 0.1.0-r1`
- `taf-rnaseq-alignment-flow = 0.1.0-r1`
- `taf-rnaseq-alignment-qc-flow = 0.1.0-r1`
- `taf-rnaseq-count-flow = 0.1.0-r1`
- `taf-rnaseq-de-flow = 0.1.0-r2`
- `taf-rnaseq-enrichment-flow = 0.1.0-r3`
- `taf-rnaseq-report-flow = 0.1.0-r3`

Each subflow records its own tool-level dependencies in its output
`versions.tsv`; the top-level flow merges those records into
`04_reports/versions.tsv`.

## Scope

This r1 release keeps Salmon as the default because it is lighter and directly
estimates transcript abundance before gene-level summarization. The optional
`--route both` branch adds splice-aware genome alignment, BAM-level QC, and
featureCounts gene counting for users who need alignments, alignment QC
evidence, or a classical alignment-count DE route.

`--route both --de-source salmon` runs the extra branch for evidence but keeps
the default Salmon DE input. `--route both --de-source featurecounts` switches
the DE input to featureCounts. The flow records the chosen route and count
source in `flow_summary.tsv`, `commands.sh`, and `run.manifest.json`.

The flow can build a Kallisto index with `--indexer both`, but r1 still uses
Salmon for expression quantification. It does not infer strandedness or
biological design automatically; set `--library-type`, `--rna-strandness`,
`--count-strand`, and `--sequencing-protocol` according to the library
preparation protocol when needed.

For example data, use `taf-rnaseq-yeast-get-data` to create a local yeast
SNF2 data directory:

```sh
taf-rnaseq-yeast-get-data \
  --outdir yeast-snf2-data-v1 \
  --stage all \
  --resume true
```

Then use `yeast-snf2-data-v1/03_results` as the data root for manual
examples. The bundled formal test can also be pointed at such a data root with
`TAFFISH_RNASEQ_TESTDATA`.
