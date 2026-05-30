# rnaseq-standard-flow De Novo Analysis Manual

This manual focuses on the explicit no-reference route of
`rnaseq-standard-flow 0.2.0-r2`. It is written for projects where a reliable
reference genome and gene annotation are missing, incomplete, or not suitable
for the biological question.

For the full standard-flow manual, including reference-mode analysis, see
`docs/manual.en.md`.

## 1. When To Use De Novo Mode

Use `--mode denovo` when the main biological object cannot be represented well
by a trusted reference genome and annotation. Common examples include
non-model organisms, mixed or poorly annotated strains, species with no usable
public genome, or projects where transcript discovery is more important than
mapping reads to known reference genes.

Do not use de novo mode only because a reference file is temporarily missing.
Missing `--genome` or `--annotation` can be a user error, so the flow never
switches to no-reference analysis automatically. You must opt in with:

```sh
--mode denovo
```

The de novo route answers a different question from the reference route:

```text
Which assembled transcripts are expressed differently, and what biological
functions can be inferred from homology evidence?
```

It does not directly answer:

```text
Which curated reference genes changed?
```

unless you later provide a reliable transcript-to-gene or transcript-cluster
mapping.

## 2. Biological and Technical Route

The no-reference standard route is:

```text
FASTQ + metadata + protein FASTA + protein-to-GO map
-> rnaseq-denovo-assembly-flow
-> rnaseq-denovo-expression-flow
-> rnaseq-denovo-annotation-flow
-> rnaseq-de-flow
-> rnaseq-enrichment-flow
-> rnaseq-report-flow
```

The biological meaning of each step is:

- Assembly reconstructs transcript sequences from reads without using a genome.
- Expression quantification estimates read support for assembled transcripts.
- ORF prediction and homology annotation connect transcript-derived proteins to known protein evidence.
- Differential expression compares transcript-level counts across biological conditions.
- Enrichment tests whether homology-derived GO gene sets are overrepresented among changed transcripts.
- The final report separates no-reference evidence from reference-mode gene-level evidence.

This route is useful, but its interpretation is more evidence-weighted than a
curated reference analysis. Assembly quality, sample depth, transcript
fragmentation, isoform redundancy, protein database choice, and GO mapping
quality all affect the final biological story.

## 3. Minimal Command

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

Add `--threads`, `--max-memory`, and project-specific design options for real
datasets:

```sh
taf-rnaseq-standard-flow \
  --mode denovo \
  --samples samples.tsv \
  --metadata metadata.tsv \
  --design '~ batch + condition' \
  --contrast condition:treated:control \
  --protein-db proteins.faa \
  --go-map protein_go_map.tsv \
  --threads 16 \
  --max-memory 64G \
  --project-name "No-reference RNA-seq" \
  --outdir rnaseq-denovo-standard-out
```

## 4. Required Inputs

### 4.1 FASTQ Sample Table

`--samples samples.tsv` is a tab-separated table with one row per biological
sample.

Required columns:

```text
sample_id	read1
```

Optional paired-end column:

```text
sample_id	read1	read2
WT_01	reads/WT_01_R1.fq.gz	reads/WT_01_R2.fq.gz
WT_02	reads/WT_02_R1.fq.gz	reads/WT_02_R2.fq.gz
KO_01	reads/KO_01_R1.fq.gz	reads/KO_01_R2.fq.gz
KO_02	reads/KO_02_R1.fq.gz	reads/KO_02_R2.fq.gz
```

Single-end and paired-end samples cannot be mixed in the same run. Relative
FASTQ paths are resolved relative to the directory containing `samples.tsv`.

For de novo assembly, balanced and biologically representative reads matter.
If one condition has much deeper sequencing, the assembler may preferentially
recover transcripts from that condition. The flow exposes Trinity read
normalization control through `--no-normalize`, but the default keeps Trinity's
normalization behavior.

### 4.2 Metadata Table

`--metadata metadata.tsv` is consumed by `rnaseq-de-flow` and must contain a
sample identifier column. By default, that column is named `sample`.

Example:

```text
sample	condition	batch
WT_01	WT	B1
WT_02	WT	B2
KO_01	KO	B1
KO_02	KO	B2
```

The sample names must match `sample_id` in `samples.tsv`. If you subset FASTQ
samples, subset metadata at the same time.

### 4.3 Design and Contrast

`--design` is the DESeq2 design formula. The simplest form is:

```sh
--design '~ condition'
```

If the experiment has a known batch variable, include it:

```sh
--design '~ batch + condition'
```

`--contrast` selects the comparison:

```sh
--contrast condition:treated:control
```

The numerator is the condition whose log2 fold-change is reported as positive
relative to the denominator.

### 4.4 Protein Database

`--protein-db proteins.faa` is a local protein FASTA used by
`rnaseq-denovo-annotation-flow` for DIAMOND homology evidence.

This is not sequencing output and is not usually delivered by a sequencing
provider. It is an external protein knowledge base prepared by the analyst for
functional annotation. For real projects, choose a high-quality protein set
from a close, well annotated species or clade, rather than a distant or overly
broad database that may turn homology annotation into noise.

Choose a database that matches the project:

- For a non-model species with a close relative, use proteins from the closest well annotated species or clade.
- For yeast-like data, a curated yeast protein set is better than a broad unrelated database.
- For broad exploratory work, a larger protein database may increase sensitivity but also runtime and ambiguous hits.

Header IDs in this FASTA must match the `subject_id` values in `--go-map`.
For example:

```text
>P12345 hypothetical kinase
MSTNPKPQR...
```

then `protein_go_map.tsv` should use `P12345` as `subject_id`.

### 4.5 Protein-to-GO Map

`--go-map protein_go_map.tsv` maps protein subject IDs to GO terms.

This is not sequencing output either. It is the functional mapping paired with
`--protein-db`; `subject_id` must match the IDs DIAMOND can hit in the protein
FASTA. The full standard-flow de novo route needs this map because enrichment
is built from the generated GMT/background resources.

Recommended columns:

```text
subject_id	go_id	go_name	namespace
P12345	GO:0004672	protein kinase activity	molecular_function
P12345	GO:0006468	protein phosphorylation	biological_process
P67890	GO:0005737	cytoplasm	cellular_component
```

Required columns:

- `subject_id`
- `go_id`

Optional but recommended:

- `go_name`
- `namespace`

The de novo annotation flow uses DIAMOND hits plus this map to generate:

```text
03_results/denovo_annotation/03_results/gene_sets/denovo_go.gmt
03_results/denovo_annotation/03_results/gene_sets/denovo_background.tsv
```

Those generated files are then used by enrichment. If the GO map has poor
coverage, enrichment may be sparse or biased.

## 5. De Novo Parameters

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--mode` | `reference` | Must be set to `denovo` for no-reference analysis. |
| `--assembler` | `trinity` | `trinity` or `rnaspades`; assembler used by `rnaseq-denovo-assembly-flow`. |
| `--max-memory` | `4G` | Memory limit passed to the assembler. Increase for real datasets. |
| `--min-contig-len` | `200` | Minimum assembled transcript length retained. |
| `--ss-lib-type` | `none` | Strand-specific library type passed to assembly where supported. |
| `--no-normalize` | off | Adds Trinity `--no_normalize_reads`; mainly for tiny or pre-normalized tests. |
| `--denovo-min-orf-aa` | `50` | Minimum TransDecoder ORF amino-acid length. |
| `--denovo-evalue` | `1e-5` | DIAMOND blastp e-value cutoff. |
| `--denovo-max-target-seqs` | `1` | DIAMOND subject hits retained per predicted protein. |

Practical choices:

- Start with `--assembler trinity` unless you have a reason to prefer rnaSPAdes.
- Increase `--max-memory` according to sample count, read depth, and organism complexity.
- Keep `--min-contig-len 200` for general use; increasing it removes short fragments but may lose small RNAs or partial transcripts.
- Use `--ss-lib-type` only when the library preparation is known.
- Keep `--denovo-max-target-seqs 1` for clean enrichment resources; increase only if you need broader annotation evidence and are prepared to handle ambiguity.

## 6. Parameters Shared With Reference Mode

`--threads` controls subflow threading. De novo assembly can be expensive, so
real projects often need more than the default `2`.

`--trim` runs fastp trimming inside de novo subflows. Use it when raw reads
still contain adapter or quality-tail problems. If reads have already been
cleaned by a trusted upstream process, leaving it off is acceptable.

`--skip-fastqc` skips FastQC. This is not recommended for first-pass project
QC, but it can be useful when upstream QC has already been archived.

`--padj-cutoff`, `--lfc-cutoff`, `--fit-type`, `--lfc-shrink`, `--min-count`,
`--min-samples`, `--top-var`, and `--top-heatmap` are passed to
`rnaseq-de-flow`. They have the same statistical meaning as in reference mode,
but the features are assembled transcripts.

`--enrichment-min-size`, `--enrichment-max-size`,
`--enrichment-pvalue-cutoff`, `--enrichment-padj-method`,
`--enrichment-top-n`, and `--enrichment-seed` are passed to
`rnaseq-enrichment-flow`. In de novo mode, the GMT/background are generated
from homology annotation.

## 7. Options Not Used In De Novo Mode

These options are reference-only:

- `--genome`
- `--annotation`
- `--gene-sets`
- `--background`
- `--route both`
- `--de-source featurecounts`
- `--rna-strandness`
- `--alignment-min-mapq`
- `--count-strand`
- `--count-feature-type`
- `--count-attribute`
- `--count-min-assigned-reads`
- `--qc-mapq`
- `--infer-sample-size`
- `--java-mem-size`
- `--sequencing-protocol`

The reason is biological as much as technical: HISAT2 alignment, featureCounts
gene counting, RSeQC, and Qualimap require a reference genome and annotation.
The de novo route instead works with assembled transcript sequences and
transcript-level matrices.

## 8. Running The Yeast Test Data In De Novo Mode

The yeast SNF2 data produced by `rnaseq-yeast-get-data` can be reused for both
the reference-mode tutorial and the de novo tutorial. First prepare the shared
data root in the same way as the reference manual:

```sh
taf-rnaseq-yeast-get-data \
  --outdir yeast-snf2-data-v1 \
  --stage all \
  --resume true

DATA="$PWD/yeast-snf2-data-v1/03_results"
```

This downloads the data once:

- The reference example uses genome, annotation, gene sets, and FASTQ files under the same `DATA` directory.
- The de novo example does not pass the genome or annotation into the flow. It reuses FASTQ files and derives a small test protein database and protein-to-GO map from the same SGD source tarball.

Yeast has an excellent reference genome, so real yeast projects should usually
prefer reference mode. This de novo example is mainly an engineering and
interpretation test for the no-reference route.

### 8.1 Fast 2v2 Subset Test

For a faster local run, select 2 WT and 2 SNF2KO samples and keep the first
50,000 reads per sample:

```sh
WORK="$PWD/yeast-denovo-test-inputs"
mkdir -p "$WORK/reads" "$WORK/resources"

{
  printf 'sample_id\tread1\tcondition\n'
  awk -F '\t' -v OFS='\t' '
    NR == 1 {
      for (i = 1; i <= NF; i++) col[$i] = i
      next
    }
    $(col["condition"]) == "snf2_KO" && ko < 2 {
      print $(col["sample_id"]), $(col["read1"]), $(col["condition"])
      ko++
      next
    }
    $(col["condition"]) == "WT" && wt < 2 {
      print $(col["sample_id"]), $(col["read1"]), $(col["condition"])
      wt++
      next
    }
  ' "$DATA/yeast-snf2-fastq-mini-v1/samples.tsv"
} > "$WORK/selected.tsv"

printf 'sample_id\tread1\tcondition\n' > "$WORK/samples.tsv"
printf 'sample\tcondition\n' > "$WORK/metadata.tsv"

tail -n +2 "$WORK/selected.tsv" |
while IFS="$(printf '\t')" read -r sample_id read1 condition; do
  case "$read1" in
    /*) fq="$read1" ;;
    *) fq="$DATA/yeast-snf2-fastq-mini-v1/$read1" ;;
  esac
  gzip -cd "$fq" \
    | awk 'NR <= 200000 { print }' > "$WORK/reads/$sample_id.fq"
  printf '%s\t%s\t%s\n' "$sample_id" "reads/$sample_id.fq" "$condition" >> "$WORK/samples.tsv"
  printf '%s\t%s\n' "$sample_id" "$condition" >> "$WORK/metadata.tsv"
done
```

Then derive the protein FASTA and protein-to-GO map needed by the de novo
annotation step from the same SGD reference source package:

```sh
REFERENCE_TAR="$DATA/yeast-reference-sgd-r64.4.1-v1/source/S288C_reference_genome_R64-4-1_20230830.tgz"
GO_TERMS="$DATA/yeast-sgd-go-gene-sets-r64.4.1-v1/metadata/go_terms.tsv"
PROTEIN_DB="$WORK/resources/yeast_orf_trans_all_R64-4-1.protein.faa"
GO_MAP="$WORK/resources/yeast_sgd_go_map.tsv"

tar -xOzf "$REFERENCE_TAR" \
  S288C_reference_genome_R64-4-1_20230830/orf_trans_all_R64-4-1_20230830.fasta.gz \
  | gzip -cd \
  | awk '
      /^>/ { print; next }
      {
        gsub(/\*/, "")
        gsub(/[[:space:]]/, "")
        if ($0 != "") print
      }
    ' > "$PROTEIN_DB"

tar -xOzf "$REFERENCE_TAR" \
  S288C_reference_genome_R64-4-1_20230830/gene_association_R64-4-1_20230830.sgd.gz \
  | gzip -cd \
  | awk -F '\t' -v OFS='\t' -v terms="$GO_TERMS" '
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
    ' > "$GO_MAP"
```

Run the de novo standard flow:

```sh
taf-rnaseq-standard-flow \
  --mode denovo \
  --samples "$WORK/samples.tsv" \
  --metadata "$WORK/metadata.tsv" \
  --design '~ condition' \
  --contrast condition:snf2_KO:WT \
  --protein-db "$PROTEIN_DB" \
  --go-map "$GO_MAP" \
  --outdir yeast-denovo-standard-out \
  --threads 4 \
  --max-memory 8G \
  --min-contig-len 100 \
  --denovo-min-orf-aa 30 \
  --denovo-evalue 1e-3 \
  --denovo-max-target-seqs 1 \
  --no-normalize \
  --project-name "Yeast SNF2 RNA-seq de novo standard test"
```

After the run, inspect:

```text
yeast-denovo-standard-out/03_results/denovo_assembly/
yeast-denovo-standard-out/03_results/denovo_expression/
yeast-denovo-standard-out/03_results/denovo_annotation/
yeast-denovo-standard-out/03_results/de/
yeast-denovo-standard-out/03_results/enrichment/
yeast-denovo-standard-out/04_reports/rnaseq_report.html
yeast-denovo-standard-out/04_reports/report_interpretation.html
```

### 8.2 Full 24-Sample Test

For the full yeast SNF2 mini package, meaning 24 biological samples with about
500,000 reads per sample, do not create a subset `WORK/samples.tsv` and do not
truncate FASTQ files. To make this section directly copy-pastable, first define
the data root, resource directory, protein DB, and GO map variables. If you
already generated `PROTEIN_DB` and `GO_MAP` in section 8.1, you can skip the
resource-generation commands and only confirm that the variables point to the
existing files.

Here, `DATA` points to the example data root produced by
`rnaseq-yeast-get-data`; the FASTQ files are the read inputs. `PROTEIN_DB` and
`GO_MAP` are annotation resources derived from the same SGD reference package
for this tutorial. They stand in for the external protein database and GO
mapping that a real no-reference project must prepare. Replace them with
species-appropriate resources for real non-model analyses.

```sh
DATA="$PWD/yeast-snf2-data-v1/03_results"
WORK="$PWD/yeast-denovo-24sample-resources"
mkdir -p "$WORK/resources"

REFERENCE_TAR="$DATA/yeast-reference-sgd-r64.4.1-v1/source/S288C_reference_genome_R64-4-1_20230830.tgz"
GO_TERMS="$DATA/yeast-sgd-go-gene-sets-r64.4.1-v1/metadata/go_terms.tsv"
PROTEIN_DB="$WORK/resources/yeast_orf_trans_all_R64-4-1.protein.faa"
GO_MAP="$WORK/resources/yeast_sgd_go_map.tsv"

tar -xOzf "$REFERENCE_TAR" \
  S288C_reference_genome_R64-4-1_20230830/orf_trans_all_R64-4-1_20230830.fasta.gz \
  | gzip -cd \
  | awk '
      /^>/ { print; next }
      {
        gsub(/\*/, "")
        gsub(/[[:space:]]/, "")
        if ($0 != "") print
      }
    ' > "$PROTEIN_DB"

tar -xOzf "$REFERENCE_TAR" \
  S288C_reference_genome_R64-4-1_20230830/gene_association_R64-4-1_20230830.sgd.gz \
  | gzip -cd \
  | awk -F '\t' -v OFS='\t' -v terms="$GO_TERMS" '
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
    ' > "$GO_MAP"
```

Then pass the package-level sample table and metadata table directly:

```sh
taf-rnaseq-standard-flow \
  --mode denovo \
  --samples "$DATA/yeast-snf2-fastq-mini-v1/samples.tsv" \
  --metadata "$DATA/yeast-snf2-fastq-mini-v1/metadata.tsv" \
  --design '~ condition' \
  --contrast condition:snf2_KO:WT \
  --protein-db "$PROTEIN_DB" \
  --go-map "$GO_MAP" \
  --outdir yeast-denovo-standard-24sample-out \
  --threads 8 \
  --max-memory 32G \
  --project-name "Yeast SNF2 RNA-seq de novo 24-sample standard"
```

This command runs the complete no-reference route on the full FASTQ package:

- `rnaseq-denovo-assembly-flow`: assemble a transcriptome from all 24 samples;
- `rnaseq-denovo-expression-flow`: quantify all 24 samples against the assembled transcripts;
- `rnaseq-denovo-annotation-flow`: attach ORF, homology, and GO evidence using `PROTEIN_DB` and `GO_MAP`;
- `rnaseq-de-flow`: run transcript-level DE for `snf2_KO` versus `WT`;
- `rnaseq-enrichment-flow`: use the de novo annotation-derived GMT/background for enrichment;
- `rnaseq-report-flow`: render the final bilingual HTML report.

The full 24-sample run is much slower and more memory-intensive than the 2v2
subset. On a server, start with at least 8 CPU threads and 32 GB memory; if
resources allow, increase `--threads` and `--max-memory` to `16` / `64G`.
For the full run, usually do not add `--no-normalize`; let Trinity's default
read normalization reduce assembly pressure. If reads have already been
cleaned by a trusted upstream process, you may add `--skip-fastqc` or leave
`--trim` off; if you want the flow to clean reads internally, add `--trim`.

## 9. Output Layout

Important directories:

```text
rnaseq-denovo-standard-out/
  00_inputs/
  01_logs/
  03_results/
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
    flow_summary.tsv
    subflows.tsv
    plot_files.tsv
    commands.sh
    versions.tsv
    methods.txt
  run.manifest.json
```

Key biological outputs:

- Assembled transcriptome FASTA from `denovo_assembly`.
- Transcript count and TPM matrices from `denovo_expression`.
- ORF, DIAMOND, and GO-derived annotation from `denovo_annotation`.
- Transcript-level DESeq2 tables and plots from `de`.
- ORA/GSEA enrichment tables and plots from `enrichment`.
- Bilingual static report from `04_reports/rnaseq_report.html`.

## 10. Reading Results

The count matrix is transcript-level. A row ID is an assembled transcript
feature, not necessarily a curated gene. Multiple assembled transcripts may
represent isoforms, paralogs, fragmented pieces, allelic variants, assembly
artifacts, or redundant assemblies.

The DE table should be interpreted as:

```text
assembled transcript X has higher or lower read support between conditions
```

not automatically as:

```text
known gene Y changed
```

The annotation table links transcripts to predicted ORFs and protein hits.
This is evidence transfer, not manual curation. Stronger confidence comes from
longer ORFs, good DIAMOND hits, consistent GO terms, and agreement across
biological replicates.

The enrichment results depend on the generated de novo background. If many
transcripts have no protein or GO support, the tested universe represents only
the annotatable subset of the assembly.

## 11. Quality Checks

After a run, inspect:

- Assembly summary tables: total transcripts, length distribution, N50-like metrics if present.
- BUSCO output from the assembly flow when available.
- Salmon quantification logs and mapping/assignment rates.
- DESeq2 sample PCA and sample correlation plots.
- The number of significant transcripts and whether one sample dominates.
- Enrichment background size and term coverage.
- `04_reports/versions.tsv` and `04_reports/commands.sh` for provenance.

Large numbers of very short transcripts, weak ORF recovery, poor protein hit
coverage, or extremely sparse GO mapping usually mean the downstream
biological interpretation should be conservative.

## 12. Troubleshooting

### The flow refuses `--route both`

`--route both` is reference-only. Use reference mode if genome alignment and
featureCounts evidence are required.

### The flow says `--protein-db` or `--go-map` is missing

De novo enrichment is built from homology evidence. Provide both files in
`--mode denovo`; the flow does not download them.

### Enrichment returns very few terms

Check whether DIAMOND hits cover enough transcripts and whether `subject_id`
values in `--go-map` match FASTA header IDs in `--protein-db`.

### DE has transcript IDs that are hard to interpret

That is expected in no-reference analysis. Use the annotation output to connect
transcript IDs to ORFs, protein hits, and GO terms. If the project needs
gene-like IDs, build a transcript clustering or transcript-to-gene map as an
additional downstream step.

### Assembly is slow or memory-limited

Increase `--max-memory`, use more `--threads`, check disk space, and consider
subsetting or normalizing reads before a pilot run. Do not use tiny test
settings for production data.

## 13. Boundaries

The de novo route does not:

- download protein databases, GO maps, BUSCO lineages, or other resources;
- create a curated gene annotation;
- guarantee one assembled transcript per biological gene;
- run reference-genome alignment or featureCounts;
- turn homology-derived annotation into manual curation;
- decide whether a no-reference analysis is biologically preferable for the project.

Use de novo mode as a reproducible no-reference baseline. For publication-grade
projects, review assembly quality, annotation coverage, and biological
plausibility before making final claims.
