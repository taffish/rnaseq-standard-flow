rnaseq-standard-flow 0.3.0-r1

Purpose:
  Run a TAFFISH RNA-seq standard analysis from FASTQ and local resources to
  expression matrices, DE tables, enrichment results, and a bilingual report.

Modes:
  --mode reference
      Default. Requires genome, annotation, gene sets, and sample metadata.
      Runs reference index, Salmon expression, DE, enrichment, and report.
      Use --route both to also run HISAT2, alignment QC, and featureCounts.

  --mode denovo
      Explicit no-reference route. Requires protein FASTA and GO map resources.
      Runs transcriptome assembly, transcript expression, homology annotation,
      DE, enrichment, and report. Missing --genome never auto-selects denovo.

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
  Chinese standard manual: https://github.com/taffish/rnaseq-standard-flow/blob/v0.3.0-r1/docs/manual.zh.md
  English standard manual: https://github.com/taffish/rnaseq-standard-flow/blob/v0.3.0-r1/docs/manual.en.md
  Chinese de novo manual: https://github.com/taffish/rnaseq-standard-flow/blob/v0.3.0-r1/docs/manual-denovo.zh.md
  English de novo manual: https://github.com/taffish/rnaseq-standard-flow/blob/v0.3.0-r1/docs/manual-denovo.en.md
  Repository and full option reference: https://github.com/taffish/rnaseq-standard-flow

Required common inputs:
  --samples PATH          TSV with sample_id, read1, optional read2.
                          Relative FASTQ paths are resolved from samples.tsv.
  --metadata PATH         TSV sample metadata for DESeq2. Default sample column: sample.
  --design TEXT           DESeq2 design formula, for example '~ condition'.
  --contrast A:B:C        DESeq2 contrast, for example condition:treated:control.
  --outdir PATH, -o PATH  Dedicated output directory. Existing output is refused
                          unless --force is used.

Reference-mode inputs:
  --genome PATH           Reference genome FASTA.
  --annotation PATH       GFF3 or GTF annotation; seq IDs must match genome FASTA.
  --gene-sets PATH        Offline GMT gene-set file in the DE feature ID space.
  --background PATH       Optional ORA background gene list.

De novo inputs:
  --protein-db PATH       External protein FASTA for DIAMOND annotation evidence.
  --go-map PATH           TSV mapping protein subject IDs to GO terms.
                          Required columns: subject_id, go_id.
                          Optional columns: go_name, namespace.

General options:
  --mode reference|denovo        Default: reference.
  --threads INT, -t INT          Default: 2.
  --route salmon|both            Reference only. Default: salmon.
  --de-source salmon|featurecounts
                                  Default: salmon. featurecounts requires
                                  --mode reference --route both.
  --trim                         Run fastp trimming where supported.
  --skip-fastqc                  Skip FastQC where supported.
  --force                        Replace outputs under an existing --outdir.
  --project-name TEXT            Report title. Default: RNA-seq project.

Reference expression options:
  --library-type TEXT            Salmon library type. Default: A.
  --indexer salmon|both          Build Salmon or Salmon+Kallisto index. Default: salmon.
  --kmer INT                     Kallisto/de novo Salmon index k-mer. Default: 31.
  --min-assigned-frags INT       Salmon minimum assigned fragments. Default: 10.
  --counts-from-abundance METHOD tximport method: no, scaledTPM,
                                  lengthScaledTPM, dtuScaledTPM. Default: no.

De novo options:
  --assembler trinity|rnaspades  Default: trinity.
  --max-memory SIZE              Assembler memory limit. Default: 4G.
  --min-contig-len INT           Minimum transcript length. Default: 200.
  --ss-lib-type none|F|R|FR|RF   Strand-specific assembly setting. Default: none.
  --no-normalize                 Add Trinity --no_normalize_reads.
  --denovo-min-orf-aa INT        TransDecoder ORF amino-acid cutoff. Default: 50.
  --denovo-evalue VALUE          DIAMOND blastp e-value cutoff. Default: 1e-5.
  --denovo-max-target-seqs INT   DIAMOND hits kept per protein. Default: 1.

Alignment/count route options:
  --rna-strandness none|F|R|FR|RF          HISAT2 strandness. Default: none.
  --alignment-min-mapq INT                 BAM MAPQ filter. Default: 0.
  --count-strand 0|1|2                     featureCounts strand mode. Default: 0.
  --count-feature-type NAME                Default: exon.
  --count-attribute NAME                   Default: gene_id.
  --count-min-assigned-reads INT           Default: 0.
  --qc-mapq INT                            Alignment QC MAPQ cutoff. Default: 30.
  --infer-sample-size INT                  RSeQC infer size. Default: 200000.
  --java-mem-size SIZE                     Qualimap Java memory. Default: 4G.
  --sequencing-protocol NAME               Qualimap protocol. Default: non-strand-specific.

Differential expression options:
  --sample-column NAME       Metadata sample column. Default: sample.
  --gene-column NAME         Reference count feature column. Default: gene_id.
                              Denovo mode uses transcript_id internally.
  --padj-cutoff FLOAT        Default: 0.05.
  --lfc-cutoff FLOAT         Absolute log2 fold-change cutoff. Default: 1.
  --fit-type TYPE            parametric, local, mean. Default: parametric.
  --lfc-shrink TYPE          none, ashr, apeglm. Default: none.
  --coef NAME                Required for --lfc-shrink apeglm.
  --min-count INT            Default: 1.
  --min-samples INT          Default: 2.
  --top-var INT              PCA feature count. Default: 500.
  --top-heatmap INT          Heatmap feature count. Default: 50.

Enrichment options:
  --enrichment-min-size INT          Default: 2.
  --enrichment-max-size INT          Default: 500.
  --enrichment-pvalue-cutoff FLOAT   Default: 1.
  --enrichment-padj-method METHOD    Default: BH.
  --enrichment-top-n INT             Default: 20.
  --enrichment-seed INT              Default: 1.

Advanced internal-tool passthrough:
  Optional expert syntax:
      @step-name: --native-tool-option value @:

  These slots default to empty. Normal runs do not need them. Full mapping,
  examples, and risks are in README and manuals.

  Index -> rnaseq-index-flow:
    @index-{agat-convert,gffread-gtf,gffread-transcripts,salmon-index,kallisto-index,hisat2-build}-step:

  Expression -> rnaseq-expression-flow:
    @expression-{fastqc-pe,fastqc-se,fastp-pe,fastp-se}-step:
    @expression-{salmon-quant-pe,salmon-quant-se,tximport,multiqc}-step:
    @expression-{salmon-quantmerge-counts,salmon-quantmerge-tpm}-step:

  Alignment -> rnaseq-alignment-flow:
    @alignment-{fastp-pe,fastp-se,hisat2-align-pe,hisat2-align-se,multiqc}-step:
    @alignment-samtools-{sort,index,quickcheck,flagstat,idxstats,mapq-filter,mapq-index}-step:

  Alignment QC -> rnaseq-alignment-qc-flow:
    @alignment-qc-samtools-{quickcheck,flagstat,idxstats}-step:
    @alignment-qc-rseqc-{bam-stat,infer-experiment,read-distribution}-step:
    @alignment-qc-{gffread-bed,qualimap-rnaseq-pe,qualimap-rnaseq-se,multiqc}-step:

  Count -> rnaseq-count-flow:
    @count-{samtools-quickcheck,featurecounts,multiqc}-step:

  De novo assembly -> rnaseq-denovo-assembly-flow:
    @denovo-assembly-{fastqc-pe,fastqc-se,fastp-pe,fastp-se}-step:
    @denovo-assembly-{trinity-assembly,rnaspades-assembly,busco,multiqc}-step:
    @denovo-assembly-seqkit-{filter,stats}-step:

  De novo expression -> rnaseq-denovo-expression-flow:
    @denovo-expression-{seqkit-transcript-stats,fastqc-pe,fastqc-se,multiqc}-step:
    @denovo-expression-{fastp-pe,fastp-se,salmon-denovo-index}-step:
    @denovo-expression-salmon-denovo-{quant-pe,quant-se,quantmerge-counts,quantmerge-tpm}-step:

  De novo annotation -> rnaseq-denovo-annotation-flow:
    @denovo-annotation-{seqkit-transcript-stats,transdecoder-longorfs,transdecoder-predict}-step:
    @denovo-annotation-diamond-{makedb,blastp}-step:

  DE -> rnaseq-de-flow:
    @de-{deseq2,deseq2-shrink,pca,de-plot}-step:

  Enrichment -> rnaseq-enrichment-flow:
    @enrichment-enrichment-{both-background,both,ora-background,ora,gsea}-step:
    @enrichment-render-{dotplot,extra-plots}-step:

  Example:
      @denovo-assembly-trinity-assembly-step: --normalize_max_read_cov 50 @:

Key outputs:
  <outdir>/03_results/
      reference/ expression/ de/ enrichment/ report/
      alignment/ alignment_qc/ count/           present with --route both
      denovo_assembly/ denovo_expression/ denovo_annotation/
      plots/png/ plots/pdf/

  <outdir>/04_reports/
      rnaseq_report.html, report_interpretation.html, commands.sh,
      versions.tsv, methods.txt, flow_summary.tsv, subflows.tsv,
      collected_files.tsv, plot_files.tsv.

  <outdir>/run.manifest.json
      Top-level provenance manifest.

Dependencies:
  taf-rnaseq-index-flow 0.2.0-r1
  taf-rnaseq-expression-flow 0.2.0-r1
  taf-rnaseq-alignment-flow 0.2.0-r1
  taf-rnaseq-alignment-qc-flow 0.2.0-r1
  taf-rnaseq-count-flow 0.2.0-r1
  taf-rnaseq-de-flow 0.2.0-r1
  taf-rnaseq-enrichment-flow 0.2.0-r1
  taf-rnaseq-report-flow 0.3.0-r2
  taf-rnaseq-denovo-assembly-flow 0.2.0-r1
  taf-rnaseq-denovo-expression-flow 0.2.0-r1
  taf-rnaseq-denovo-annotation-flow 0.2.0-r1

Boundaries:
  Offline umbrella flow. It does not download references, gene sets, protein
  databases, GO mappings, models, or genome resources. It does not infer
  design, strandedness, or whether reference or de novo analysis is appropriate.
  De novo transcript IDs are assembled transcript features, not known reference genes.

Wrapper options:
  -h, --help       Show this help.
  -v, --version    Show package and command version.
  --compile        Print generated shell code instead of running it.
