#!/bin/sh

set -eu

script_dir=$(CDPATH= cd "$(dirname "$0")" && pwd)
project_dir=$(CDPATH= cd "$script_dir/.." && pwd)
rnaseq_root=$(CDPATH= cd "$project_dir/.." && pwd)
default_data_root=$(CDPATH= cd "$rnaseq_root/test-data/yeast/data/03_results" 2>/dev/null && pwd || printf '%s\n' "$rnaseq_root/test-data/yeast/data/03_results")

DATA="${TAFFISH_RNASEQ_TESTDATA:-$default_data_root}"
OUT="${TAFFISH_RNASEQ_REAL_REFERENCE_OUT:-$script_dir/test-real-reference-out}"
THREADS="${TAFFISH_RNASEQ_REAL_THREADS:-32}"

DATA=$(CDPATH= cd "$DATA" 2>/dev/null && pwd || printf '%s\n' "$DATA")

TAFFISH_CONTAINER_BACKEND="${TAFFISH_CONTAINER_BACKEND:-apptainer}"
TAF_HISTORY_MODE="${TAF_HISTORY_MODE:-off}"
export TAFFISH_CONTAINER_BACKEND TAF_HISTORY_MODE

flow_cmd="$project_dir/target/taf-rnaseq-standard-flow-v0.3.0-r1"

need_file() {
    if [ ! -s "$1" ]; then
        echo "test-real-reference: missing required file: $1" >&2
        exit 1
    fi
}

samples="$DATA/yeast-snf2-fastq-mini-v1/samples.tsv"
genome="$DATA/yeast-reference-sgd-r64.4.1-v1/reference/genome/yeast_s288c_reference_genome_R64-4-1.fa"
annotation="$DATA/yeast-reference-sgd-r64.4.1-v1/reference/annotation/yeast_s288c_gene_annotation_R64-4-1.gff3"
metadata="$DATA/yeast-snf2-fastq-mini-v1/metadata.tsv"
gene_sets="$DATA/yeast-sgd-go-gene-sets-r64.4.1-v1/gene_sets/sgd_go_bp.gmt"
background="$DATA/yeast-sgd-go-gene-sets-r64.4.1-v1/background/yeast_background_genes.tsv"

need_file "$samples"
need_file "$genome"
need_file "$annotation"
need_file "$metadata"
need_file "$gene_sets"
need_file "$background"

cd "$project_dir"

echo "[REAL reference] taf check"
taf check

echo "[REAL reference] taf build"
taf build

if [ ! -x "$flow_cmd" ]; then
    echo "test-real-reference: built wrapper not executable: $flow_cmd" >&2
    exit 1
fi

echo "[REAL reference] run rnaseq-standard-flow"
"$flow_cmd" \
    --samples "$samples" \
    --genome "$genome" \
    --annotation "$annotation" \
    --metadata "$metadata" \
    --design '~ condition' \
    --contrast condition:snf2_KO:WT \
    --gene-sets "$gene_sets" \
    --background "$background" \
    --outdir "$OUT" \
    --route both \
    --de-source featurecounts \
    --threads "$THREADS" \
    --count-strand 0 \
    --count-feature-type exon \
    --count-attribute gene_id \
    --qc-mapq 30 \
    --infer-sample-size 200000 \
    --java-mem-size 128G \
    --sequencing-protocol non-strand-specific \
    --fit-type parametric \
    --min-count 10 \
    --min-samples 4 \
    --padj-cutoff 0.05 \
    --lfc-cutoff 1 \
    --top-var 500 \
    --top-heatmap 50 \
    --enrichment-min-size 2 \
    --enrichment-max-size 500 \
    --enrichment-top-n 20 \
    --project-name "Yeast SNF2 RNA-seq reference standard" \
    --force

echo "[REAL reference] done"
echo "  report: $OUT/04_reports/rnaseq_report.html"
echo "  interpretation: $OUT/04_reports/report_interpretation.html"
