#!/bin/sh
set -eu

script_dir=$(CDPATH= cd "$(dirname "$0")" && pwd)
project_dir=$(CDPATH= cd "$script_dir/.." && pwd)
bio_apps_dir=$(CDPATH= cd "$project_dir/../../.." && pwd)
rnaseq_root=$(CDPATH= cd "$project_dir/.." && pwd)

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

if ! command -v taf >/dev/null 2>&1; then
    echo "smoke: taf command not found in PATH." >&2
    exit 127
fi

if ! command -v taffish >/dev/null 2>&1; then
    echo "smoke: taffish command not found in PATH." >&2
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
    taf-rnaseq-report-flow-v0.1.0-r3
do
    if ! command -v "$dep" >/dev/null 2>&1; then
        echo "smoke: dependency wrapper not found in PATH: $dep" >&2
        exit 127
    fi
done

TAFFISH_CONTAINER_BACKEND=${TAFFISH_CONTAINER_BACKEND:-podman}
export TAFFISH_CONTAINER_BACKEND
TAF_HISTORY_MODE=${TAF_HISTORY_MODE:-off}
export TAF_HISTORY_MODE

tmpdir=$(mktemp -d "$project_dir/.taf-smoke.XXXXXX")
cleanup() {
    cd "$project_dir" 2>/dev/null || :
    rm -rf "$tmpdir"
}
trap cleanup EXIT INT TERM HUP

cd "$project_dir"

echo "[SMOKE] taf check"
taf check

echo "[SMOKE] taf build"
taf build

flow_cmd="$project_dir/target/taf-rnaseq-standard-flow-v0.1.0-r1"
if [ ! -x "$flow_cmd" ]; then
    echo "smoke: built flow command is missing or not executable: $flow_cmd" >&2
    exit 1
fi

echo "[SMOKE] help and version"
"$flow_cmd" --help >/dev/null
"$flow_cmd" --version >/dev/null

run_dir="$tmpdir/run"
fixture="$run_dir/fixture"
reads_dir="$fixture/reads"
mkdir -p "$reads_dir"

echo "[SMOKE] create synthetic RNA-seq fixture"
awk -v genome="$fixture/genome.fa" -v gff="$fixture/annotation.gff3" -v seqs="$fixture/gene_sequences.tsv" '
    function seq_for_gene(i, len,    p, x, s, bases) {
        bases = "ACGT"
        x = i * 7919 + 17
        s = ""
        for (p = 1; p <= len; p++) {
            x = (x * 1103515245 + 12345) % 2147483647
            s = s substr(bases, (x % 4) + 1, 1)
        }
        return s
    }
    BEGIN {
        print ">chr1" > genome
        print "##gff-version 3" > gff
        print "gene_id\tsequence" > seqs
        chr = ""
        for (i = 1; i <= 40; i++) {
            gene = sprintf("gene%02d", i)
            tx = sprintf("tx%02d", i)
            exon = sprintf("exon%02d", i)
            start = (i - 1) * 220 + 1
            end = start + 199
            seq = seq_for_gene(i, 200)
            chr = chr seq "NNNNNNNNNNNNNNNNNNNN"
            printf "chr1\tsmoke\tgene\t%d\t%d\t.\t+\t.\tID=%s;Name=%s\n", start, end, gene, gene > gff
            printf "chr1\tsmoke\tmRNA\t%d\t%d\t.\t+\t.\tID=%s;Parent=%s;Name=%s\n", start, end, tx, gene, tx > gff
            printf "chr1\tsmoke\texon\t%d\t%d\t.\t+\t.\tID=%s;Parent=%s\n", start, end, exon, tx > gff
            print gene "\t" seq > seqs
        }
        for (i = 1; i <= length(chr); i += 80) {
            print substr(chr, i, 80) > genome
        }
    }
'

cat > "$fixture/counts_template.tsv" <<'EOF'
gene_id	c1	c2	c3	t1	t2	t3
gene01	7	15	12	7	4	7
gene02	12	7	8	7	6	15
gene03	51	20	12	16	31	26
gene04	36	46	28	21	41	26
gene05	69	71	80	109	158	132
gene06	90	41	33	50	18	37
gene07	221	141	173	300	374	203
gene08	282	295	132	99	164	135
gene09	10	14	5	12	9	10
gene10	20	8	13	12	13	12
gene11	18	13	17	16	16	13
gene12	29	33	8	24	30	26
gene13	32	30	31	268	115	151
gene14	94	93	103	32	49	53
gene15	117	145	78	142	186	128
gene16	485	319	387	255	94	81
gene17	6	7	6	5	7	5
gene18	15	21	15	16	4	12
gene19	33	32	33	25	35	32
gene20	38	44	49	33	44	41
gene21	44	56	65	83	177	189
gene22	93	66	105	45	46	34
gene23	93	174	62	293	347	142
gene24	301	202	516	211	208	41
gene25	5	11	9	12	6	14
gene26	13	11	23	16	6	15
gene27	14	32	14	18	39	33
gene28	26	80	14	32	18	39
gene29	89	61	51	184	114	132
gene30	91	74	66	52	32	29
gene31	89	72	207	463	163	282
gene32	205	168	184	151	124	131
gene33	9	5	12	8	12	7
gene34	16	13	6	12	5	18
gene35	25	33	19	28	14	15
gene36	44	24	24	16	15	29
gene37	51	41	79	169	151	123
gene38	90	93	69	55	35	30
gene39	159	83	98	122	291	448
gene40	233	249	271	92	92	108
EOF

awk -F '\t' \
    -v samples="$fixture/samples.tsv" \
    -v metadata="$fixture/metadata.tsv" \
    -v reads_dir="$reads_dir" '
    function qual(n,    i, q) {
        q = ""
        for (i = 1; i <= n; i++) q = q "I"
        return q
    }
    FNR == NR {
        if (FNR > 1) seq[$1] = $2
        next
    }
    FNR == 1 {
        print "sample_id\tread1\tcondition" > samples
        print "sample\tcondition" > metadata
        for (i = 2; i <= NF; i++) {
            sid[i] = $i
            group[i] = ($i ~ /^c/) ? "control" : "treated"
            fq[i] = reads_dir "/" sid[i] ".fq"
            print sid[i] "\treads/" sid[i] ".fq\t" group[i] > samples
            print sid[i] "\t" group[i] > metadata
        }
        next
    }
    {
        gene = $1
        for (i = 2; i <= NF; i++) {
            n = int($i)
            for (r = 1; r <= n; r++) {
                start = 1 + ((r * 7 + i * 3) % 120)
                read = substr(seq[gene], start, 75)
                printf "@%s_%s_%03d\n%s\n+\n%s\n", sid[i], gene, r, read, qual(length(read)) > fq[i]
            }
        }
    }
' "$fixture/gene_sequences.tsv" "$fixture/counts_template.tsv"

{
    printf 'set_up\tup genes\tgene01\tgene02\tgene03\tgene04\tgene05\n'
    printf 'set_down\tdown genes\tgene06\tgene07\tgene08\tgene09\tgene10\n'
    printf 'set_background\tall genes'
    i=1
    while [ "$i" -le 40 ]; do
        printf '\tgene%02d' "$i"
        i=$((i + 1))
    done
    printf '\n'
} > "$fixture/gene_sets.gmt"

{
    printf 'gene_id\n'
    i=1
    while [ "$i" -le 40 ]; do
        printf 'gene%02d\n' "$i"
        i=$((i + 1))
    done
} > "$fixture/background.tsv"

echo "[SMOKE] rnaseq-standard-flow synthetic full route"
(
    cd "$run_dir"
    "$flow_cmd" \
        --samples "$fixture/samples.tsv" \
        --genome "$fixture/genome.fa" \
        --annotation "$fixture/annotation.gff3" \
        --metadata "$fixture/metadata.tsv" \
        --design '~ condition' \
        --contrast condition:treated:control \
        --gene-sets "$fixture/gene_sets.gmt" \
        --background "$fixture/background.tsv" \
        --outdir standard-out \
        --threads 1 \
        --skip-fastqc \
        --min-assigned-frags 1 \
        --fit-type mean \
        --padj-cutoff 1 \
        --lfc-cutoff 0 \
        --min-count 1 \
        --min-samples 2 \
        --top-var 20 \
        --top-heatmap 10 \
        --enrichment-min-size 1 \
        --enrichment-max-size 50 \
        --enrichment-top-n 5 \
        --project-name "TAFFISH standard smoke"
)
cd "$project_dir"

out="$run_dir/standard-out"

echo "[SMOKE] output checks"
test -s "$out/00_inputs/standard_inputs.tsv"
test -s "$out/01_logs/flow.log"
test -s "$out/01_logs/steps/02_rnaseq_index.log"
test -s "$out/01_logs/steps/03_rnaseq_expression.log"
test -s "$out/01_logs/steps/04_rnaseq_alignment.log"
test -s "$out/01_logs/steps/05_rnaseq_alignment_qc.log"
test -s "$out/01_logs/steps/06_rnaseq_count.log"
test -s "$out/01_logs/steps/07_rnaseq_de.log"
test -s "$out/01_logs/steps/08_rnaseq_enrichment.log"
test -s "$out/01_logs/steps/09_rnaseq_report.log"
test -s "$out/01_logs/steps/10_collect_reports.log"
test -s "$out/03_results/reference/03_results/salmon_index/info.json"
test -s "$out/03_results/reference/03_results/tx2gene.tsv"
test -s "$out/03_results/expression/03_results/matrices/gene_counts.tsv"
test -s "$out/03_results/de/03_results/de/results.tsv"
test -s "$out/03_results/de/03_results/gene_lists/significant_genes.tsv"
test -s "$out/03_results/de/03_results/gene_lists/ranked_genes.tsv"
test -s "$out/03_results/de/03_results/plots/plot_summary.tsv"
test -s "$out/03_results/enrichment/03_results/enrichment/ora_results.tsv"
test -s "$out/03_results/enrichment/03_results/enrichment/gsea_results.tsv"
test -s "$out/03_results/enrichment/03_results/enrichment/dotplot.png"
test -s "$out/03_results/enrichment/03_results/enrichment/plot_summary.tsv"
test -s "$out/03_results/enrichment/03_results/enrichment/ora_barplot.png"
test -s "$out/03_results/enrichment/03_results/enrichment/gsea_nes_plot.png"
test -s "$out/03_results/enrichment/03_results/enrichment/gsea_enrichment_curves.png"
test -s "$out/03_results/report/04_reports/rnaseq_report.html"
test -d "$out/03_results/plots"
test -s "$out/04_reports/rnaseq_report.html"
test -s "$out/04_reports/commands.sh"
test -s "$out/04_reports/versions.tsv"
test -s "$out/04_reports/methods.txt"
test -s "$out/04_reports/flow_summary.tsv"
test -s "$out/04_reports/subflows.tsv"
test -s "$out/04_reports/plot_files.tsv"
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

grep -F 'route	salmon' "$out/04_reports/flow_summary.tsv" >/dev/null
grep -F 'de_source	salmon' "$out/04_reports/flow_summary.tsv" >/dev/null
grep -F 'sample_count	6' "$out/04_reports/flow_summary.tsv" >/dev/null
grep -F 'standard_plots	28' "$out/04_reports/flow_summary.tsv" >/dev/null
grep -F 'de	deg_counts_barplot	pdf' "$out/04_reports/plot_files.tsv" >/dev/null
grep -F 'de	sample_correlation_heatmap	png' "$out/04_reports/plot_files.tsv" >/dev/null
grep -F 'enrichment	dotplot	png' "$out/04_reports/plot_files.tsv" >/dev/null
grep -F 'enrichment	ora_barplot	pdf' "$out/04_reports/plot_files.tsv" >/dev/null
grep -F 'enrichment	gsea_nes_plot	png' "$out/04_reports/plot_files.tsv" >/dev/null
grep -F 'enrichment	gsea_enrichment_curves	pdf' "$out/04_reports/plot_files.tsv" >/dev/null
grep -F 'rnaseq-index-flow' "$out/04_reports/subflows.tsv" >/dev/null
grep -F 'rnaseq-expression-flow' "$out/04_reports/subflows.tsv" >/dev/null
grep -F '04_rnaseq_alignment	rnaseq-alignment-flow		skipped' "$out/04_reports/subflows.tsv" >/dev/null
grep -F '05_rnaseq_alignment_qc	rnaseq-alignment-qc-flow		skipped' "$out/04_reports/subflows.tsv" >/dev/null
grep -F '06_rnaseq_count	rnaseq-count-flow		skipped' "$out/04_reports/subflows.tsv" >/dev/null
grep -F 'rnaseq-de-flow' "$out/04_reports/subflows.tsv" >/dev/null
grep -F 'rnaseq-enrichment-flow' "$out/04_reports/subflows.tsv" >/dev/null
grep -F 'rnaseq-report-flow' "$out/04_reports/subflows.tsv" >/dev/null
grep -F 'taf-rnaseq-index-flow-v0.1.0-r1' "$out/04_reports/commands.sh" >/dev/null
grep -F 'taf-rnaseq-expression-flow-v0.1.0-r1' "$out/04_reports/commands.sh" >/dev/null
grep -F 'taf-rnaseq-de-flow-v0.1.0-r2' "$out/04_reports/commands.sh" >/dev/null
grep -F 'taf-rnaseq-enrichment-flow-v0.1.0-r3' "$out/04_reports/commands.sh" >/dev/null
grep -F 'taf-rnaseq-report-flow-v0.1.0-r3 --standard-out' "$out/04_reports/commands.sh" >/dev/null
grep -F 'TAFFISH standard smoke' "$out/04_reports/rnaseq_report.html" >/dev/null
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
grep -F '../03_results/report/03_results/collected_plots/de.pca_plot.png' "$out/04_reports/rnaseq_report.html" >/dev/null
grep -F '../03_results/report/03_results/collected_html/expression.multiqc/index.html' "$out/04_reports/rnaseq_report.html" >/dev/null
if grep -F '../03_results/collected_' "$out/04_reports/rnaseq_report.html" >/dev/null; then
    echo "smoke: top-level report contains stale report-flow relative links." >&2
    exit 1
fi
grep -F 'provided_modules	5' "$out/03_results/report/04_reports/project_summary.tsv" >/dev/null
grep -F 'plot_groups	14' "$out/03_results/report/04_reports/project_summary.tsv" >/dev/null
grep -F '"flow": "rnaseq-standard-flow"' "$out/run.manifest.json" >/dev/null
grep -F '"route": "salmon"' "$out/run.manifest.json" >/dev/null
grep -F '"de_source": "salmon"' "$out/run.manifest.json" >/dev/null
if command -v python3 >/dev/null 2>&1; then
    python3 -m json.tool "$out/run.manifest.json" >/dev/null
fi

echo "[SMOKE] existing outdir is refused"
if (
    cd "$run_dir"
    "$flow_cmd" \
        --samples "$fixture/samples.tsv" \
        --genome "$fixture/genome.fa" \
        --annotation "$fixture/annotation.gff3" \
        --metadata "$fixture/metadata.tsv" \
        --design '~ condition' \
        --contrast condition:treated:control \
        --gene-sets "$fixture/gene_sets.gmt" \
        --background "$fixture/background.tsv" \
        --outdir standard-out \
        --skip-fastqc
) >/dev/null 2>&1; then
    echo "smoke: existing outdir was not refused." >&2
    exit 1
fi

echo "[SMOKE] --force rerun"
(
    cd "$run_dir"
    "$flow_cmd" \
        --samples "$fixture/samples.tsv" \
        --genome "$fixture/genome.fa" \
        --annotation "$fixture/annotation.gff3" \
        --metadata "$fixture/metadata.tsv" \
        --design '~ condition' \
        --contrast condition:treated:control \
        --gene-sets "$fixture/gene_sets.gmt" \
        --background "$fixture/background.tsv" \
        --outdir standard-out \
        --route both \
        --de-source featurecounts \
        --threads 1 \
        --skip-fastqc \
        --min-assigned-frags 1 \
        --alignment-min-mapq 0 \
        --count-strand 0 \
        --count-feature-type exon \
        --count-attribute gene_id \
        --count-min-assigned-reads 0 \
        --qc-mapq 0 \
        --infer-sample-size 1000 \
        --java-mem-size 2G \
        --sequencing-protocol non-strand-specific \
        --fit-type mean \
        --padj-cutoff 1 \
        --lfc-cutoff 0 \
        --min-count 1 \
        --min-samples 2 \
        --enrichment-min-size 1 \
        --enrichment-max-size 50 \
        --force
)
test -s "$out/04_reports/rnaseq_report.html"
test -s "$out/03_results/reference/03_results/hisat2_index/genome.1.ht2"
test -s "$out/03_results/alignment/04_reports/bam_files.tsv"
test -s "$out/03_results/alignment_qc/04_reports/rnaseq_qc_summary.tsv"
test -s "$out/03_results/count/03_results/matrices/gene_counts.tsv"
test -s "$out/03_results/plots/png/de.sample_correlation_heatmap.png"
test -s "$out/03_results/plots/png/enrichment.dotplot.png"
test -s "$out/03_results/plots/png/enrichment.gsea_nes_plot.png"
test -s "$out/03_results/plots/pdf/de.sample_correlation_heatmap.pdf"
test -s "$out/03_results/plots/pdf/enrichment.dotplot.pdf"
test -s "$out/03_results/plots/pdf/enrichment.gsea_enrichment_curves.pdf"
test -s "$out/04_reports/plot_files.tsv"
grep -F 'route	both' "$out/04_reports/flow_summary.tsv" >/dev/null
grep -F 'de_source	featurecounts' "$out/04_reports/flow_summary.tsv" >/dev/null
grep -F 'standard_plots	28' "$out/04_reports/flow_summary.tsv" >/dev/null
grep -F '04_rnaseq_alignment	rnaseq-alignment-flow' "$out/04_reports/subflows.tsv" >/dev/null
grep -F '05_rnaseq_alignment_qc	rnaseq-alignment-qc-flow' "$out/04_reports/subflows.tsv" >/dev/null
grep -F '06_rnaseq_count	rnaseq-count-flow' "$out/04_reports/subflows.tsv" >/dev/null
grep -F 'taf-rnaseq-alignment-flow-v0.1.0-r1' "$out/04_reports/commands.sh" >/dev/null
grep -F 'taf-rnaseq-alignment-qc-flow-v0.1.0-r1' "$out/04_reports/commands.sh" >/dev/null
grep -F 'taf-rnaseq-count-flow-v0.1.0-r1' "$out/04_reports/commands.sh" >/dev/null
grep -F 'taf-rnaseq-report-flow-v0.1.0-r3 --standard-out' "$out/04_reports/commands.sh" >/dev/null
grep -F 'provided_modules	8' "$out/03_results/report/04_reports/project_summary.tsv" >/dev/null
grep -F 'plot_groups	14' "$out/03_results/report/04_reports/project_summary.tsv" >/dev/null
grep -F '"route": "both"' "$out/run.manifest.json" >/dev/null
grep -F '"de_source": "featurecounts"' "$out/run.manifest.json" >/dev/null
if command -v python3 >/dev/null 2>&1; then
    python3 -m json.tool "$out/run.manifest.json" >/dev/null
fi

stray=$(find "$run_dir" -mindepth 1 -maxdepth 1 ! -name fixture ! -name standard-out -print)
if [ -n "$stray" ]; then
    echo "smoke: flow wrote unexpected files outside outdir:" >&2
    printf '%s\n' "$stray" >&2
    exit 1
fi

echo "[SMOKE] ok"
