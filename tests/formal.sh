#!/bin/sh
set -eu

script_dir=$(CDPATH= cd "$(dirname "$0")" && pwd)
project_dir=$(CDPATH= cd "$script_dir/.." && pwd)
bio_apps_dir=$(CDPATH= cd "$project_dir/../../.." && pwd)
rnaseq_root=$(CDPATH= cd "$project_dir/.." && pwd)
default_data_root=$(CDPATH= cd "$rnaseq_root/test-data/yeast/data/03_results" 2>/dev/null && pwd || printf '%s\n' "$rnaseq_root/test-data/yeast/data/03_results")
data_root=${TAFFISH_RNASEQ_TESTDATA:-$default_data_root}

for target_dir in \
    "$rnaseq_root/subflows/rnaseq-index-flow/target" \
    "$rnaseq_root/subflows/rnaseq-expression-flow/target" \
    "$rnaseq_root/subflows/rnaseq-alignment-flow/target" \
    "$rnaseq_root/subflows/rnaseq-alignment-qc-flow/target" \
    "$rnaseq_root/subflows/rnaseq-count-flow/target" \
    "$rnaseq_root/subflows/rnaseq-de-flow/target" \
    "$rnaseq_root/subflows/rnaseq-enrichment-flow/target" \
    "$rnaseq_root/subflows/rnaseq-report-flow/target" \
    "$bio_apps_dir/tools/agat/target" \
    "$bio_apps_dir/tools/gffread/target" \
    "$bio_apps_dir/tools/salmon/target" \
    "$bio_apps_dir/tools/kallisto/target" \
    "$bio_apps_dir/tools/hisat2/target" \
    "$bio_apps_dir/tools/samtools/target" \
    "$bio_apps_dir/tools/rseqc/target" \
    "$bio_apps_dir/tools/qualimap/target" \
    "$bio_apps_dir/tools/subread/target" \
    "$bio_apps_dir/tools/fastqc/target" \
    "$bio_apps_dir/tools/fastp/target" \
    "$bio_apps_dir/tools/multiqc/target" \
    "$bio_apps_dir/tools/bioconductor-rnaseq/target" \
    "$bio_apps_dir/tools/enrichment-r/target"
do
    if [ -d "$target_dir" ]; then
        PATH="$target_dir:$PATH"
    fi
done
export PATH

TAFFISH_CONTAINER_BACKEND=${TAFFISH_CONTAINER_BACKEND:-podman}
export TAFFISH_CONTAINER_BACKEND
TAF_HISTORY_MODE=${TAF_HISTORY_MODE:-off}
export TAF_HISTORY_MODE

skip_formal() {
    echo "formal: skipped: $*" >&2
    exit 0
}

if [ ! -d "$data_root" ]; then
    skip_formal "RNA-seq formal data root not found: $data_root"
fi

fastq_samples="$data_root/yeast-snf2-fastq-mini-v1/samples.tsv"
genome="$data_root/yeast-reference-sgd-r64.4.1-v1/reference/genome/yeast_s288c_reference_genome_R64-4-1.fa"
annotation="$data_root/yeast-reference-sgd-r64.4.1-v1/reference/annotation/yeast_s288c_gene_annotation_R64-4-1.gff3"
gene_sets_pkg="$data_root/yeast-sgd-go-gene-sets-r64.4.1-v1"
gene_sets="$gene_sets_pkg/gene_sets/sgd_go_bp.gmt"
background="$gene_sets_pkg/background/yeast_background_genes.tsv"

[ -s "$fastq_samples" ] || skip_formal "missing yeast FASTQ samples.tsv: $fastq_samples"
[ -s "$genome" ] || skip_formal "missing yeast reference genome FASTA: $genome"
[ -s "$annotation" ] || skip_formal "missing yeast reference GFF3 annotation: $annotation"
[ -s "$gene_sets" ] || skip_formal "missing yeast SGD GO BP GMT: $gene_sets"
[ -s "$background" ] || skip_formal "missing yeast enrichment background: $background"

if ! command -v taf >/dev/null 2>&1; then
    echo "formal: taf command not found in PATH." >&2
    exit 127
fi

for dep in \
    taf-rnaseq-index-flow-v0.1.0-r1 \
    taf-rnaseq-expression-flow-v0.1.0-r1 \
    taf-rnaseq-alignment-flow-v0.1.0-r1 \
    taf-rnaseq-alignment-qc-flow-v0.1.0-r1 \
    taf-rnaseq-count-flow-v0.1.0-r1 \
    taf-rnaseq-de-flow-v0.1.0-r2 \
    taf-rnaseq-enrichment-flow-v0.1.0-r3 \
    taf-rnaseq-report-flow-v0.1.0-r4
do
    if ! command -v "$dep" >/dev/null 2>&1; then
        echo "formal: dependency wrapper not found in PATH: $dep" >&2
        exit 127
    fi
done

tmpdir=$(mktemp -d "$project_dir/.taf-formal.XXXXXX")
cleanup() {
    cd "$project_dir" 2>/dev/null || :
    rm -rf "$tmpdir"
}
trap cleanup EXIT INT TERM HUP

cd "$project_dir"

echo "[FORMAL] taf check"
taf check

echo "[FORMAL] taf build"
taf build

flow_cmd="$project_dir/target/taf-rnaseq-standard-flow-v0.1.0-r2"
if [ ! -x "$flow_cmd" ]; then
    echo "formal: built flow command is missing or not executable: $flow_cmd" >&2
    exit 1
fi

run_dir="$tmpdir/run"
mkdir -p "$run_dir"

formal_samples="$run_dir/samples.subset.tsv"
formal_metadata="$run_dir/metadata.subset.tsv"
awk -F '\t' -v OFS='\t' -v base="$(dirname "$fastq_samples")" -v metadata="$formal_metadata" '
    NR == 1 {
        for (i = 1; i <= NF; i++) col[$i] = i
        if (!("sample_id" in col) || !("read1" in col) || !("condition" in col)) {
            print "formal: samples table must contain sample_id, read1, and condition" > "/dev/stderr"
            exit 2
        }
        print "sample_id", "read1", "condition"
        print "sample", "condition" > metadata
        next
    }
    $0 == "" { next }
    $(col["condition"]) == "snf2_KO" && ko < 2 {
        print $(col["sample_id"]), base "/" $(col["read1"]), $(col["condition"])
        print $(col["sample_id"]), $(col["condition"]) > metadata
        ko++
        next
    }
    $(col["condition"]) == "WT" && wt < 2 {
        print $(col["sample_id"]), base "/" $(col["read1"]), $(col["condition"])
        print $(col["sample_id"]), $(col["condition"]) > metadata
        wt++
        next
    }
    END {
        close(metadata)
        if (ko < 2 || wt < 2) {
            print "formal: need at least 2 snf2_KO and 2 WT samples" > "/dev/stderr"
            exit 3
        }
    }
' "$fastq_samples" > "$formal_samples"

echo "[FORMAL] rnaseq-standard-flow yeast FASTQ subset"
(
    cd "$run_dir"
    "$flow_cmd" \
        --samples "$formal_samples" \
        --genome "$genome" \
        --annotation "$annotation" \
        --metadata "$formal_metadata" \
        --design '~ condition' \
        --contrast condition:snf2_KO:WT \
        --gene-sets "$gene_sets" \
        --background "$background" \
        --outdir standard-out \
        --route both \
        --de-source featurecounts \
        --threads 2 \
        --skip-fastqc \
        --min-assigned-frags 1 \
        --alignment-min-mapq 0 \
        --count-strand 0 \
        --count-feature-type exon \
        --count-attribute gene_id \
        --count-min-assigned-reads 0 \
        --qc-mapq 30 \
        --infer-sample-size 200000 \
        --java-mem-size 4G \
        --sequencing-protocol non-strand-specific \
        --fit-type mean \
        --min-count 1 \
        --min-samples 2 \
        --padj-cutoff 1 \
        --lfc-cutoff 0 \
        --top-var 100 \
        --top-heatmap 30 \
        --enrichment-min-size 2 \
        --enrichment-max-size 500 \
        --enrichment-top-n 20 \
        --project-name "Yeast SNF2 RNA-seq standard formal"
)

out="$run_dir/standard-out"
test -s "$out/03_results/reference/03_results/salmon_index/info.json"
test -s "$out/03_results/reference/03_results/tx2gene.tsv"
test -s "$out/03_results/reference/03_results/hisat2_index/genome.1.ht2"
test -s "$out/03_results/expression/03_results/matrices/gene_counts.tsv"
test -s "$out/03_results/expression/03_results/matrices/gene_tpm.tsv"
test -s "$out/03_results/alignment/04_reports/bam_files.tsv"
test -s "$out/03_results/alignment_qc/04_reports/rnaseq_qc_summary.tsv"
test -s "$out/03_results/count/03_results/matrices/gene_counts.tsv"
test -s "$out/03_results/de/03_results/de/results.tsv"
test -s "$out/03_results/de/03_results/gene_lists/significant_genes.tsv"
test -s "$out/03_results/de/03_results/gene_lists/ranked_genes.tsv"
test -s "$out/03_results/de/03_results/plots/plot_summary.tsv"
test -s "$out/03_results/enrichment/03_results/enrichment/ora_results.tsv"
test -s "$out/03_results/enrichment/03_results/enrichment/gsea_results.tsv"
test -s "$out/03_results/enrichment/03_results/enrichment/dotplot.png"
test -s "$out/03_results/enrichment/03_results/enrichment/dotplot.original.png"
test -s "$out/03_results/enrichment/03_results/enrichment/plot_summary.tsv"
test -s "$out/03_results/enrichment/03_results/enrichment/ora_barplot.png"
test -s "$out/03_results/enrichment/03_results/enrichment/gsea_nes_plot.png"
test -s "$out/03_results/enrichment/03_results/enrichment/gsea_enrichment_curves.png"
test -s "$out/03_results/report/04_reports/rnaseq_report.html"
test -d "$out/03_results/plots"
test -s "$out/04_reports/rnaseq_report.html"
test -s "$out/04_reports/project_summary.tsv"
test -s "$out/04_reports/collected_files.tsv"
test -s "$out/04_reports/plot_files.tsv"
test -s "$out/04_reports/commands.sh"
test -s "$out/04_reports/versions.tsv"
test -s "$out/04_reports/methods.txt"
test -s "$out/04_reports/flow_summary.tsv"
test -s "$out/04_reports/subflows.tsv"
test -s "$out/run.manifest.json"

for plot_file in \
    pdf/de.pca_plot.pdf \
    png/de.pca_plot.png \
    pdf/de.ma_plot.pdf \
    png/de.ma_plot.png \
    pdf/de.volcano_plot.pdf \
    png/de.volcano_plot.png \
    pdf/de.deg_counts_barplot.pdf \
    png/de.deg_counts_barplot.png \
    pdf/de.heatmap.pdf \
    png/de.heatmap.png \
    pdf/de.sample_correlation_heatmap.pdf \
    png/de.sample_correlation_heatmap.png \
    pdf/de.expression_distribution.pdf \
    png/de.expression_distribution.png \
    pdf/de.normalized_count_distribution.pdf \
    png/de.normalized_count_distribution.png \
    pdf/de.top_genes_expression.pdf \
    png/de.top_genes_expression.png \
    pdf/enrichment.dotplot.pdf \
    png/enrichment.dotplot.png \
    pdf/enrichment.dotplot_original.pdf \
    png/enrichment.dotplot_original.png \
    pdf/enrichment.ora_barplot.pdf \
    png/enrichment.ora_barplot.png \
    pdf/enrichment.gsea_nes_plot.pdf \
    png/enrichment.gsea_nes_plot.png \
    pdf/enrichment.gsea_enrichment_curves.pdf \
    png/enrichment.gsea_enrichment_curves.png
do
    test -s "$out/03_results/plots/$plot_file"
done

grep -F 'sample_count	4' "$out/04_reports/flow_summary.tsv" >/dev/null
grep -F 'route	both' "$out/04_reports/flow_summary.tsv" >/dev/null
grep -F 'de_source	featurecounts' "$out/04_reports/flow_summary.tsv" >/dev/null
grep -F 'standard_plots	28' "$out/04_reports/flow_summary.tsv" >/dev/null
grep -F 'de	top_genes_expression	png' "$out/04_reports/plot_files.tsv" >/dev/null
grep -F 'enrichment	dotplot_original	pdf' "$out/04_reports/plot_files.tsv" >/dev/null
grep -F 'enrichment	ora_barplot	png' "$out/04_reports/plot_files.tsv" >/dev/null
grep -F 'enrichment	gsea_nes_plot	pdf' "$out/04_reports/plot_files.tsv" >/dev/null
grep -F 'enrichment	gsea_enrichment_curves	png' "$out/04_reports/plot_files.tsv" >/dev/null
grep -F 'Yeast SNF2 RNA-seq standard formal' "$out/04_reports/rnaseq_report.html" >/dev/null
grep -F 'taf-rnaseq-index-flow-v0.1.0-r1' "$out/04_reports/commands.sh" >/dev/null
grep -F 'taf-rnaseq-expression-flow-v0.1.0-r1' "$out/04_reports/commands.sh" >/dev/null
grep -F 'taf-rnaseq-alignment-flow-v0.1.0-r1' "$out/04_reports/commands.sh" >/dev/null
grep -F 'taf-rnaseq-alignment-qc-flow-v0.1.0-r1' "$out/04_reports/commands.sh" >/dev/null
grep -F 'taf-rnaseq-count-flow-v0.1.0-r1' "$out/04_reports/commands.sh" >/dev/null
grep -F 'taf-rnaseq-de-flow-v0.1.0-r2' "$out/04_reports/commands.sh" >/dev/null
grep -F 'taf-rnaseq-enrichment-flow-v0.1.0-r3' "$out/04_reports/commands.sh" >/dev/null
grep -F 'taf-rnaseq-report-flow-v0.1.0-r4 --standard-out' "$out/04_reports/commands.sh" >/dev/null
grep -F 'TAFFISH RNA-seq project report' "$out/04_reports/rnaseq_report.html" >/dev/null
grep -F 'data-lang-toggle="zh"' "$out/04_reports/rnaseq_report.html" >/dev/null
grep -F 'https://taffish.github.io/' "$out/04_reports/rnaseq_report.html" >/dev/null
grep -F 'IntersectionObserver' "$out/04_reports/rnaseq_report.html" >/dev/null
grep -F 'workflow-map' "$out/04_reports/rnaseq_report.html" >/dev/null
grep -F 'Deliverables and Output Structure' "$out/04_reports/rnaseq_report.html" >/dev/null
grep -F 'Readable enrichment dotplot' "$out/04_reports/rnaseq_report.html" >/dev/null
grep -F 'ORA top-term barplot' "$out/04_reports/rnaseq_report.html" >/dev/null
grep -F 'GSEA NES ranking' "$out/04_reports/rnaseq_report.html" >/dev/null
grep -F 'GSEA enrichment curves' "$out/04_reports/rnaseq_report.html" >/dev/null
grep -F 'Top gene expression' "$out/04_reports/rnaseq_report.html" >/dev/null
test -s "$out/04_reports/report_interpretation.html"
grep -F 'guide-sidebar' "$out/04_reports/report_interpretation.html" >/dev/null
grep -F 'Long-Form Module Chapters' "$out/04_reports/report_interpretation.html" >/dev/null
grep -F '长文模块章节' "$out/04_reports/report_interpretation.html" >/dev/null
grep -F '../03_results/report/03_results/collected_plots/de.pca_plot.png' "$out/04_reports/rnaseq_report.html" >/dev/null
grep -F '../03_results/report/03_results/collected_html/expression.multiqc/index.html' "$out/04_reports/rnaseq_report.html" >/dev/null
if grep -F '../03_results/collected_' "$out/04_reports/rnaseq_report.html" >/dev/null; then
    echo "formal: top-level report contains stale report-flow relative links." >&2
    exit 1
fi
grep -F 'provided_modules	8' "$out/03_results/report/04_reports/project_summary.tsv" >/dev/null
grep -F 'plot_groups	14' "$out/03_results/report/04_reports/project_summary.tsv" >/dev/null
grep -F 'rnaseq-index-flow' "$out/04_reports/subflows.tsv" >/dev/null
grep -F 'rnaseq-expression-flow' "$out/04_reports/subflows.tsv" >/dev/null
grep -F 'rnaseq-alignment-flow' "$out/04_reports/subflows.tsv" >/dev/null
grep -F 'rnaseq-alignment-qc-flow' "$out/04_reports/subflows.tsv" >/dev/null
grep -F 'rnaseq-count-flow' "$out/04_reports/subflows.tsv" >/dev/null
grep -F 'rnaseq-de-flow' "$out/04_reports/subflows.tsv" >/dev/null
grep -F 'rnaseq-enrichment-flow' "$out/04_reports/subflows.tsv" >/dev/null
grep -F 'rnaseq-report-flow' "$out/04_reports/subflows.tsv" >/dev/null
grep -F '"flow": "rnaseq-standard-flow"' "$out/run.manifest.json" >/dev/null
grep -F '"route": "both"' "$out/run.manifest.json" >/dev/null
grep -F '"de_source": "featurecounts"' "$out/run.manifest.json" >/dev/null

if command -v python3 >/dev/null 2>&1; then
    python3 -m json.tool "$out/run.manifest.json" >/dev/null
fi

echo "[FORMAL] ok"
