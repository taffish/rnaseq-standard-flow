#!/bin/sh

set -eu

export TAFFISH_CONTAINER_BACKEND=podman
export TAF_HISTORY_MODE=off

DATA="${TAFFISH_RNASEQ_TESTDATA:-../../test-data/yeast/data/03_results}"
OUT="${TAFFISH_RNASEQ_REAL_OUT:-./test-real-run-out}"
THREADS="${TAFFISH_RNASEQ_REAL_THREADS:-6}"

taf build

../target/taf-rnaseq-standard-flow-v0.1.0-r1 \
  --samples "$DATA/yeast-snf2-fastq-mini-v1/samples.tsv" \
  --genome "$DATA/yeast-reference-sgd-r64.4.1-v1/reference/genome/yeast_s288c_reference_genome_R64-4-1.fa" \
  --annotation "$DATA/yeast-reference-sgd-r64.4.1-v1/reference/annotation/yeast_s288c_gene_annotation_R64-4-1.gff3" \
  --metadata "$DATA/yeast-snf2-fastq-mini-v1/metadata.tsv" \
  --design '~ condition' \
  --contrast condition:snf2_KO:WT \
  --gene-sets "$DATA/yeast-sgd-go-gene-sets-r64.4.1-v1/gene_sets/sgd_go_bp.gmt" \
  --background "$DATA/yeast-sgd-go-gene-sets-r64.4.1-v1/background/yeast_background_genes.tsv" \
  --outdir "$OUT" \
  --route both \
  --threads "$THREADS" \
  --count-strand 0 \
  --count-feature-type exon \
  --count-attribute gene_id \
  --qc-mapq 30 \
  --infer-sample-size 200000 \
  --java-mem-size 4G \
  --sequencing-protocol non-strand-specific \
  --fit-type local \
  --min-count 10 \
  --min-samples 4 \
  --padj-cutoff 0.05 \
  --lfc-cutoff 1 \
  --top-var 500 \
  --top-heatmap 50 \
  --enrichment-min-size 2 \
  --enrichment-max-size 500 \
  --enrichment-top-n 20 \
  --project-name "Yeast SNF2 RNA-seq standard manual"
