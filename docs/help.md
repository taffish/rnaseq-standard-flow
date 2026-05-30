rnaseq-standard-flow 0.2.0-r2

Purpose:
  Run a TAFFISH RNA-seq standard analysis from FASTQ and local resources to
  expression matrices, differential expression, enrichment, and a bilingual
  static project report.

Compatibility:
  0.2.0-r2 keeps the previous reference-mode interface and behavior unchanged by default.
  Existing reference commands that pass --samples, --genome, --annotation,
  --metadata, --design, --contrast, --gene-sets, and --outdir continue to run
  as the reference route.

Modes:
  --mode reference
      Default. Use genome + annotation, build reference indexes, quantify with
      Salmon, then run DE, enrichment, and report. Optional --route both adds
      HISAT2 alignment, alignment QC, and featureCounts counting.

  --mode denovo
      Explicit no-reference route. Assemble a transcriptome, quantify reads
      against assembled transcripts, transfer homology/GO evidence from local
      protein resources, then run transcript-level DE, enrichment, and report.
      This mode must be selected explicitly; missing --genome does not trigger
      automatic de novo analysis.

Reference usage:
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

De novo usage:
  taf-rnaseq-standard-flow \
    --mode denovo \
    --samples samples.tsv \
    --metadata metadata.tsv \
    --design '~ condition' \
    --contrast condition:treated:control \
    --protein-db proteins.faa \
    --go-map protein_go_map.tsv \
    --outdir rnaseq-denovo-standard-out

Detailed manuals:
  Chinese standard manual: https://github.com/taffish/rnaseq-standard-flow/blob/v0.2.0-r2/docs/manual.zh.md
  English standard manual: https://github.com/taffish/rnaseq-standard-flow/blob/v0.2.0-r2/docs/manual.en.md
  Chinese de novo manual: https://github.com/taffish/rnaseq-standard-flow/blob/v0.2.0-r2/docs/manual-denovo.zh.md
  English de novo manual: https://github.com/taffish/rnaseq-standard-flow/blob/v0.2.0-r2/docs/manual-denovo.en.md
  Repository: https://github.com/taffish/rnaseq-standard-flow

Required inputs for all modes:
  --samples PATH
      Tab-separated FASTQ sample table. Required columns are sample_id and
      read1. Optional read2 enables paired reads. Relative FASTQ paths are
      resolved relative to samples.tsv.

  --metadata PATH
      Tab-separated sample metadata consumed by rnaseq-de-flow. By default it
      must contain a sample column named sample; change with --sample-column.

  --design TEXT
      DESeq2 design formula, for example: '~ condition'.

  --contrast FACTOR:NUMERATOR:DENOMINATOR
      DESeq2 contrast, for example condition:treated:control.

  --outdir PATH, -o PATH
      Dedicated output directory. Existing directories are refused unless
      --force is set.

Reference-mode required inputs:
  --genome PATH
      Reference genome FASTA consumed by rnaseq-index-flow.

  --annotation PATH
      GFF3 or GTF annotation consumed by rnaseq-index-flow. Sequence IDs must
      match the first token of genome FASTA headers.

  --gene-sets PATH
      Offline GMT gene-set file consumed by rnaseq-enrichment-flow.

Reference-mode optional input:
  --background PATH
      Optional background gene list for ORA enrichment. When provided, the
      flow filters DE gene-list and ranked-gene inputs to this ID space before
      calling rnaseq-enrichment-flow and records
      04_reports/enrichment_background_filter.tsv.

De novo required inputs:
  --protein-db PATH
      Local protein FASTA database for rnaseq-denovo-annotation-flow DIAMOND
      evidence. No database is downloaded.
      This is not sequencing output. It is an external annotation resource
      chosen by the analyst, usually from a close annotated species or clade.

  --go-map PATH
      TSV mapping protein subject IDs to GO terms. Header columns are
      subject_id, go_id, and optional go_name, namespace. The de novo route
      uses this to generate denovo_go.gmt and denovo_background.tsv.
      This is not sequencing output either. It must match --protein-db:
      subject_id values should correspond to protein FASTA IDs hit by DIAMOND.

Main options:
  --mode reference|denovo
      Analysis mode. Default: reference.

  --threads INT
      Threads for subflows. Default: 2.

  --route salmon|both
      Reference-mode analysis route. Default: salmon.
      both also builds a HISAT2 index, aligns reads, runs alignment QC, and
      counts genes with featureCounts. In 0.2.0-r2, --mode denovo supports only
      --route salmon.

  --de-source salmon|featurecounts
      Count matrix used by rnaseq-de-flow. Default: salmon.
      featurecounts requires --mode reference --route both. In denovo mode,
      DE uses transcript_counts.tsv and --gene-column transcript_id.

  --trim
      Run fastp trimming inside expression or de novo subflows.

  --skip-fastqc
      Skip FastQC inside expression or de novo subflows.

Reference expression options:
  --library-type TEXT
      Salmon library type. Default: A.

  --indexer salmon|both
      Reference index mode passed to rnaseq-index-flow. Default: salmon.
      both also builds a Kallisto index; reference expression quantification still uses Salmon.

  --min-assigned-frags INT
      Salmon minimum assigned fragments. Default: 10.

  --counts-from-abundance METHOD
      tximport gene-count handling: no, scaledTPM, lengthScaledTPM, or
      dtuScaledTPM. Default: no.

De novo options:
  --assembler trinity|rnaspades
      Assembler for rnaseq-denovo-assembly-flow. Default: trinity.

  --max-memory SIZE
      Assembler memory limit, for example 4G or 32G. Default: 4G.

  --min-contig-len INT
      Minimum assembled transcript length retained. Default: 200.

  --ss-lib-type none|F|R|FR|RF
      Strand-specific library type passed to assembly where supported.
      Default: none.

  --no-normalize
      Add Trinity --no_normalize_reads through rnaseq-denovo-assembly-flow.
      Useful for tiny tests and pre-downsampled datasets.

  --denovo-min-orf-aa INT
      Minimum ORF amino-acid length for TransDecoder. Default: 50.

  --denovo-evalue VALUE
      DIAMOND blastp e-value threshold. Default: 1e-5.

  --denovo-max-target-seqs INT
      Maximum DIAMOND target sequences kept per predicted protein. Default: 1.

Alignment/count route options:
  --rna-strandness none|F|R|FR|RF
      HISAT2 RNA strandness. Default: none.

  --alignment-min-mapq INT
      MAPQ filter for alignment and counting. Default: 0.

  --count-strand 0|1|2
      featureCounts strand mode. Default: 0.

  --count-feature-type NAME
      featureCounts feature type. Default: exon.

  --count-attribute NAME
      featureCounts grouping attribute. Default: gene_id.

  --count-min-assigned-reads INT
      Minimum assigned reads required by rnaseq-count-flow. Default: 0.

  --qc-mapq INT
      MAPQ cutoff for alignment QC. Default: 30.

  --infer-sample-size INT
      RSeQC infer_experiment.py sample size. Default: 200000.

  --java-mem-size SIZE
      Qualimap Java memory setting. Default: 4G.

  --sequencing-protocol NAME
      Qualimap protocol: non-strand-specific, strand-specific-forward, or
      strand-specific-reverse. Default: non-strand-specific.

Differential expression options:
  --sample-column NAME         Metadata sample column. Default: sample.
  --gene-column NAME           Reference-mode count matrix ID column.
                                Default: gene_id. Denovo mode fixes this to
                                transcript_id.
  --padj-cutoff FLOAT          Adjusted-P cutoff for significant features.
                                Default: 0.05.
  --lfc-cutoff FLOAT           Absolute log2 fold-change cutoff. Default: 1.
  --fit-type TYPE              DESeq2 dispersion fit: parametric, local, mean.
                                Default: parametric.
  --lfc-shrink TYPE            none, ashr, or apeglm. Default: none.
  --coef NAME                  Required only for --lfc-shrink apeglm.
  --min-count INT              Minimum count filter. Default: 1.
  --min-samples INT            Minimum sample filter. Default: 2.
  --top-var INT                Features used for PCA selection. Default: 500.
  --top-heatmap INT            Features shown in heatmap. Default: 50.

Enrichment options:
  --enrichment-min-size INT        Minimum GMT set size. Default: 2.
  --enrichment-max-size INT        Maximum GMT set size. Default: 500.
  --enrichment-pvalue-cutoff FLOAT Enrichment P-value cutoff. Default: 1.
  --enrichment-padj-method METHOD  Adjustment method. Default: BH.
  --enrichment-top-n INT           Top sets used for plots. Default: 20.
  --enrichment-seed INT            Seed passed to enrichment-r. Default: 1.

Report options:
  --project-name TEXT
      Project name shown in the static report. Default: RNA-seq project.

  --force
      Replace rnaseq-standard-flow outputs under an existing --outdir.

Outputs:
  <outdir>/00_inputs/
      Input table/resource snapshots and standard_inputs.tsv.

  <outdir>/01_logs/
      Flow log and per-step logs for each subflow call.

  <outdir>/03_results/reference/
  <outdir>/03_results/expression/
      Reference-mode rnaseq-index-flow and rnaseq-expression-flow outputs.

  <outdir>/03_results/alignment/
  <outdir>/03_results/alignment_qc/
  <outdir>/03_results/count/
      Optional reference-mode alignment/count branch outputs when
      --route both is used.

  <outdir>/03_results/denovo_assembly/
  <outdir>/03_results/denovo_expression/
  <outdir>/03_results/denovo_annotation/
      De novo mode outputs: assembled transcriptome, transcript count/TPM
      matrices, homology annotation, denovo_go.gmt, and denovo background.

  <outdir>/03_results/de/
      rnaseq-de-flow output. Reference mode uses gene-level counts by default;
      denovo mode uses transcript-level counts.

  <outdir>/03_results/enrichment/
      rnaseq-enrichment-flow ORA/GSEA tables and plots.

  <outdir>/03_results/report/
      Full rnaseq-report-flow 0.2.0-r2 collector output.

  <outdir>/03_results/plots/png/
  <outdir>/03_results/plots/pdf/
      Standard-flow DE and enrichment plot collection split by format.

  <outdir>/04_reports/
      Top-level bilingual rnaseq_report.html, report_interpretation.html,
      commands.sh, versions.tsv, methods.txt, flow_summary.tsv, subflows.tsv,
      collected_files.tsv, plot_files.tsv, and optional
      enrichment_background_filter.tsv. The report receives the selected mode,
      route, and DE count source explicitly, so overview labels distinguish
      reference and de novo analysis correctly.

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
  taf-rnaseq-report-flow 0.2.0-r2
  taf-rnaseq-denovo-assembly-flow 0.1.0-r1
  taf-rnaseq-denovo-expression-flow 0.1.0-r1
  taf-rnaseq-denovo-annotation-flow 0.1.0-r1

Boundaries:
  0.2.0-r2 is an offline umbrella flow. It does not download references, gene sets,
  protein databases, GO mappings, models, or genome resources at runtime. It
  does not infer biological design, strandedness, or whether a project should
  be reference-based or de novo. De novo transcript IDs are assembled
  transcript features, not known reference genes, and homology-derived
  annotation is evidence rather than manually curated function.

Wrapper options:
  -h, --help       Show this help.
  -v, --version    Show package and command version.
  --compile        Print generated shell code instead of running it.
