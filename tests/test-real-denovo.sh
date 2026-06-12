#!/bin/sh

set -eu

script_dir=$(CDPATH= cd "$(dirname "$0")" && pwd)
project_dir=$(CDPATH= cd "$script_dir/.." && pwd)
rnaseq_root=$(CDPATH= cd "$project_dir/.." && pwd)
default_data_root=$(CDPATH= cd "$rnaseq_root/test-data/yeast/data/03_results" 2>/dev/null && pwd || printf '%s\n' "$rnaseq_root/test-data/yeast/data/03_results")

DATA="${TAFFISH_RNASEQ_TESTDATA:-$default_data_root}"
OUT="${TAFFISH_RNASEQ_REAL_DENOVO_OUT:-$script_dir/test-real-denovo-out}"
WORK="${TAFFISH_RNASEQ_REAL_DENOVO_WORK:-$script_dir/test-real-denovo-inputs}"
THREADS="${TAFFISH_RNASEQ_REAL_THREADS:-64}"
MAX_MEMORY="${TAFFISH_RNASEQ_REAL_DENOVO_MAX_MEMORY:-256G}"
TRINITY_NORM_COV="${TAFFISH_RNASEQ_REAL_DENOVO_TRINITY_NORM_COV:-50}"

DATA=$(CDPATH= cd "$DATA" 2>/dev/null && pwd || printf '%s\n' "$DATA")

TAFFISH_CONTAINER_BACKEND="${TAFFISH_CONTAINER_BACKEND:-apptainer}"
TAF_HISTORY_MODE="${TAF_HISTORY_MODE:-off}"
export TAFFISH_CONTAINER_BACKEND TAF_HISTORY_MODE

flow_cmd="$project_dir/target/taf-rnaseq-standard-flow-v0.3.0-r1"

need_file() {
    if [ ! -s "$1" ]; then
        echo "test-real-denovo: missing required file: $1" >&2
        exit 1
    fi
}

samples="$DATA/yeast-snf2-fastq-mini-v1/samples.tsv"
metadata="$DATA/yeast-snf2-fastq-mini-v1/metadata.tsv"
reference_tar="$DATA/yeast-reference-sgd-r64.4.1-v1/source/S288C_reference_genome_R64-4-1_20230830.tgz"
go_terms="$DATA/yeast-sgd-go-gene-sets-r64.4.1-v1/metadata/go_terms.tsv"

need_file "$samples"
need_file "$metadata"
need_file "$reference_tar"
need_file "$go_terms"

mkdir -p "$WORK/resources"
protein_db="$WORK/resources/yeast_orf_trans_all_R64-4-1.protein.faa"
go_map="$WORK/resources/yeast_sgd_go_map.tsv"

echo "[REAL denovo] prepare protein DB"
tar -xOzf "$reference_tar" S288C_reference_genome_R64-4-1_20230830/orf_trans_all_R64-4-1_20230830.fasta.gz \
    | gzip -cd \
    | awk '
        /^>/ { print; next }
        {
            gsub(/\*/, "")
            gsub(/[[:space:]]/, "")
            if ($0 != "") print
        }
    ' > "$protein_db"
need_file "$protein_db"

echo "[REAL denovo] prepare protein-to-GO map"
tar -xOzf "$reference_tar" S288C_reference_genome_R64-4-1_20230830/gene_association_R64-4-1_20230830.sgd.gz \
    | gzip -cd \
    | awk -F '\t' -v OFS='\t' -v terms="$go_terms" '
        BEGIN {
            while ((getline line < terms) > 0) {
                n = split(line, t, "\t")
                if (n < 3 || t[1] == "go_id") continue
                go_name[t[1]] = t[2]
                go_namespace[t[1]] = t[3]
            }
            close(terms)
            print "subject_id", "go_id", "go_name", "namespace"
        }
        /^!/ || NF < 11 { next }
        {
            split($11, synonyms, "|")
            subject = synonyms[1]
            go = $5
            if (subject == "" || go == "") next
            name = (go in go_name) ? go_name[go] : go
            namespace = (go in go_namespace) ? go_namespace[go] : $9
            if (namespace == "P") namespace = "biological_process"
            else if (namespace == "F") namespace = "molecular_function"
            else if (namespace == "C") namespace = "cellular_component"
            key = subject SUBSEP go
            if (!(key in seen)) {
                seen[key] = 1
                print subject, go, name, namespace
            }
        }
    ' > "$go_map"
need_file "$go_map"

cd "$project_dir"

echo "[REAL denovo] taf check"
taf check

echo "[REAL denovo] taf build"
taf build

if [ ! -x "$flow_cmd" ]; then
    echo "test-real-denovo: built wrapper not executable: $flow_cmd" >&2
    exit 1
fi

echo "[REAL denovo] run rnaseq-standard-flow"
"$flow_cmd" \
    --mode denovo \
    --samples "$samples" \
    --metadata "$metadata" \
    --design '~ condition' \
    --contrast condition:snf2_KO:WT \
    --protein-db "$protein_db" \
    --go-map "$go_map" \
    --outdir "$OUT" \
    --threads "$THREADS" \
    --max-memory "$MAX_MEMORY" \
    --min-contig-len 200 \
    --denovo-min-orf-aa 50 \
    --denovo-evalue 1e-5 \
    --denovo-max-target-seqs 1 \
    --project-name "Yeast SNF2 RNA-seq de novo standard" \
    @denovo-assembly-trinity-assembly-step: --normalize_max_read_cov "$TRINITY_NORM_COV" @: \
    --force

echo "[REAL denovo] done"
echo "  report: $OUT/04_reports/rnaseq_report.html"
echo "  interpretation: $OUT/04_reports/report_interpretation.html"
