rnaseq-standard-flow 0.1.0-r2

Purpose:
  Run the TAFFISH RNA-seq r2 standard route from FASTQ/reference inputs to a
  static project report. The default r2 route is Salmon-first:

    rnaseq-index-flow -> rnaseq-expression-flow -> rnaseq-de-flow
    -> rnaseq-enrichment-flow -> rnaseq-report-flow

  With --route both, the flow also runs the alignment/count branch:

    rnaseq-index-flow --genome-indexer hisat2
    -> rnaseq-alignment-flow -> rnaseq-alignment-qc-flow
    -> rnaseq-count-flow

  Differential expression uses Salmon/tximport counts by default. Set
  --route both --de-source featurecounts to run DESeq2 from the featureCounts
  matrix instead.

Detailed manuals:
  Chinese: https://github.com/taffish/rnaseq-standard-flow/blob/v0.1.0-r2/docs/manual.zh.md
  English: https://github.com/taffish/rnaseq-standard-flow/blob/v0.1.0-r2/docs/manual.en.md
  Repository: https://github.com/taffish/rnaseq-standard-flow

Usage:
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

Required inputs:
  --samples PATH
      Tab-separated FASTQ sample table consumed by rnaseq-expression-flow.
      Required columns are sample_id and read1. Optional read2 enables paired
      reads. Relative FASTQ paths are resolved relative to this samples.tsv.

  --genome PATH
      Reference genome FASTA consumed by rnaseq-index-flow.

  --annotation PATH
      GFF3 or GTF annotation consumed by rnaseq-index-flow. Sequence IDs must
      match the first token of the genome FASTA headers.

  --metadata PATH
      Tab-separated sample metadata consumed by rnaseq-de-flow. By default it
      must contain a sample column named sample; change with --sample-column.

  --design TEXT
      DESeq2 design formula, for example: '~ condition'.

  --contrast FACTOR:NUMERATOR:DENOMINATOR
      DESeq2 contrast, for example condition:treated:control.

  --gene-sets PATH
      Offline GMT gene-set file consumed by rnaseq-enrichment-flow.

  --outdir PATH, -o PATH
      Dedicated output directory. Existing directories are refused unless
      --force is set.

Optional inputs:
  --background PATH
      Optional background gene list for ORA enrichment. When provided,
      rnaseq-standard-flow filters DE gene-list and ranked-gene inputs to this
      background ID space before calling rnaseq-enrichment-flow, and records
      the filtering summary in 04_reports/enrichment_background_filter.tsv.

Main options:
  --threads INT
      Threads for index and expression subflows. Default: 2.

  --route salmon|both
      Analysis route. Default: salmon.
      salmon runs reference, Salmon expression, DE, enrichment, and report.
      both also builds a HISAT2 index, aligns reads, runs alignment QC, and
      counts genes with featureCounts.

  --de-source salmon|featurecounts
      Count matrix used by rnaseq-de-flow. Default: salmon.
      featurecounts requires --route both.

  --library-type TEXT
      Salmon library type passed to rnaseq-expression-flow. Default: A.

  --indexer salmon|both
      Reference index mode passed to rnaseq-index-flow. Default: salmon.
      The r2 standard route consumes the Salmon index. both also builds the
      Kallisto index as part of the reference bundle.

  --trim
      Run fastp trimming inside rnaseq-expression-flow before Salmon.

  --skip-fastqc
      Skip FastQC inside rnaseq-expression-flow.

  --min-assigned-frags INT
      Salmon minimum assigned fragments. Default: 10.

  --counts-from-abundance METHOD
      tximport gene-count handling: no, scaledTPM, lengthScaledTPM, or
      dtuScaledTPM. Default: no.

Alignment/count route options:
  --rna-strandness none|F|R|FR|RF
      HISAT2 RNA strandness passed to rnaseq-alignment-flow. Default: none.

  --alignment-min-mapq INT
      MAPQ filter passed to rnaseq-alignment-flow and rnaseq-count-flow.
      Default: 0.

  --count-strand 0|1|2
      featureCounts strand mode passed to rnaseq-count-flow. Default: 0.

  --count-feature-type NAME
      featureCounts feature type, usually exon. Default: exon.

  --count-attribute NAME
      featureCounts grouping attribute, usually gene_id. Default: gene_id.

  --count-min-assigned-reads INT
      Minimum assigned reads required by rnaseq-count-flow. Default: 0.

  --qc-mapq INT
      MAPQ cutoff for RSeQC in rnaseq-alignment-qc-flow. Default: 30.

  --infer-sample-size INT
      RSeQC infer_experiment.py sample size. Default: 200000.

  --java-mem-size SIZE
      Qualimap Java memory setting. Default: 4G.

  --sequencing-protocol NAME
      Qualimap protocol: non-strand-specific, strand-specific-forward, or
      strand-specific-reverse. Default: non-strand-specific.

Differential expression options:
  --sample-column NAME         Metadata sample column. Default: sample.
  --gene-column NAME           Count matrix gene column. Default: gene_id.
  --padj-cutoff FLOAT          Adjusted-P cutoff for significant genes.
                                Default: 0.05.
  --lfc-cutoff FLOAT           Absolute log2 fold-change cutoff. Default: 1.
  --fit-type TYPE              DESeq2 dispersion fit: parametric, local, mean.
                                Default: parametric.
  --lfc-shrink TYPE            none, ashr, or apeglm. Default: none.
  --coef NAME                  Required only for --lfc-shrink apeglm.
  --min-count INT              Minimum count filter. Default: 1.
  --min-samples INT            Minimum sample filter. Default: 2.
  --top-var INT                Genes used for PCA selection. Default: 500.
  --top-heatmap INT            Genes shown in heatmap. Default: 50.

Enrichment options:
  --enrichment-min-size INT
      Minimum GMT set size. Default: 2.

  --enrichment-max-size INT
      Maximum GMT set size. Default: 500.

  --enrichment-pvalue-cutoff FLOAT
      Enrichment P-value cutoff. Default: 1.

  --enrichment-padj-method METHOD
      Adjustment method: holm, hochberg, hommel, bonferroni, BH, BY, fdr, or
      none. Default: BH.

  --enrichment-top-n INT
      Number of top sets used for enrichment plots. Default: 20.

  --enrichment-seed INT
      Seed passed to enrichment-r. Default: 1.

Report options:
  --project-name TEXT
      Project name shown in the static report. Default: RNA-seq project.

  --force
      Replace rnaseq-standard-flow outputs under an existing --outdir.

Outputs:
  <outdir>/00_inputs/
      Snapshots of small input tables and standard_inputs.tsv.

  <outdir>/01_logs/
      Flow log and per-step logs for each subflow call.

  <outdir>/03_results/reference/
      rnaseq-index-flow output, including Salmon index and tx2gene.tsv.

  <outdir>/03_results/expression/
      rnaseq-expression-flow output, including Salmon quant files and gene
      expression matrices.

  <outdir>/03_results/alignment/
      Optional rnaseq-alignment-flow output when --route both is used,
      including sorted BAM files and bam_files.tsv.

  <outdir>/03_results/alignment_qc/
      Optional rnaseq-alignment-qc-flow output when --route both is used,
      including SAMtools/RSeQC/Qualimap summaries.

  <outdir>/03_results/count/
      Optional rnaseq-count-flow output when --route both is used, including
      featureCounts results and gene_counts.tsv.

  <outdir>/03_results/de/
      rnaseq-de-flow output, including DESeq2 tables, plots, and gene lists.
      It uses Salmon counts by default, or featureCounts counts when
      --route both --de-source featurecounts is set.

  <outdir>/03_results/enrichment/
      rnaseq-enrichment-flow output, including ORA/GSEA tables and plots.

  <outdir>/03_results/report/
      Full rnaseq-report-flow r4 collector output.

  <outdir>/03_results/plots/
      Standard-flow plot collection. This includes DE r2 and enrichment r3
      plots in both PDF and PNG form, split into png/ and pdf/ subdirectories.
      Figures include PCA, MA, volcano, DEG counts, heatmap, sample
      correlation, expression distributions, top genes expression, and
      enrichment dotplots, ORA barplot, GSEA NES plot, and GSEA enrichment
      curves.

  <outdir>/04_reports/
      Top-level r4 bilingual rnaseq_report.html and report_interpretation.html,
      commands.sh, versions.tsv, methods.txt, flow_summary.tsv, subflows.tsv,
      collected_files.tsv, plot_files.tsv, and optional
      enrichment_background_filter.tsv. The HTML report includes workflow
      diagrams, active sidebar navigation, deliverables/output-structure notes,
      linked QC/report bundles, and organized DE/ORA/GSEA plot sections. The
      interpretation guide adds a floating contents sidebar and long-form
      RNA-seq biology/technology explanations.

  <outdir>/run.manifest.json
      Standard-flow provenance manifest.

Dependencies:
  taf-rnaseq-index-flow 0.1.0-r1
  taf-rnaseq-expression-flow 0.1.0-r1
  taf-rnaseq-alignment-flow 0.1.0-r1
  taf-rnaseq-alignment-qc-flow 0.1.0-r1
  taf-rnaseq-count-flow 0.1.0-r1
  taf-rnaseq-de-flow 0.1.0-r2
  taf-rnaseq-enrichment-flow 0.1.0-r3
  taf-rnaseq-report-flow 0.1.0-r4

Boundaries:
  r2 is an offline umbrella flow with a default Salmon route and an explicit
  optional HISAT2/featureCounts route. It does not download references, gene
  sets, genome databases, or annotation resources at runtime. It does not infer
  biological design, choose strandedness automatically, remove failed samples,
  or replace a full production workflow engine.
