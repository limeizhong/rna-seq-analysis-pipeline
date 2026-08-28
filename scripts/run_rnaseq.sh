#!/usr/bin/env bash
# Conventional paired-end RNA-seq: fastp -> HISAT2 -> featureCounts -> DESeq2.
set -euo pipefail

usage() { echo "Usage: $0 --config FILE --mode full|counts" >&2; exit 2; }
CONFIG=""; MODE=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --config) CONFIG=${2:-}; shift 2 ;;
    --mode) MODE=${2:-}; shift 2 ;;
    -h|--help) usage ;;
    *) echo "Unknown option: $1" >&2; usage ;;
  esac
done
[[ -n "$CONFIG" && -f "$CONFIG" && ( "$MODE" == full || "$MODE" == counts ) ]] || usage
# shellcheck disable=SC1090
source "$CONFIG"
SKILL_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_NAME=${PROJECT_NAME:-RNA-seq项目}; SPECIES_NAME=${SPECIES_NAME:-未说明物种}; REFERENCE_VERSION=${REFERENCE_VERSION:-未说明版本}
WARN_MIN_REPLICATES=${WARN_MIN_REPLICATES:-3}; WARN_MIN_READ_RETENTION=${WARN_MIN_READ_RETENTION:-70}; WARN_MIN_Q30=${WARN_MIN_Q30:-80}; WARN_MIN_ALIGNMENT_RATE=${WARN_MIN_ALIGNMENT_RATE:-70}; WARN_MIN_ASSIGNMENT_RATE=${WARN_MIN_ASSIGNMENT_RATE:-50}; WARN_MIN_FILTERED_GENE_FRACTION=${WARN_MIN_FILTERED_GENE_FRACTION:-0.05}; WARN_MAX_GROUP_LIBRARY_DEPTH_RATIO=${WARN_MAX_GROUP_LIBRARY_DEPTH_RATIO:-2}; WARN_MIN_WITHIN_GROUP_CORRELATION=${WARN_MIN_WITHIN_GROUP_CORRELATION:-0.80}; WARN_MIN_ENRICHMENT_MAPPING_RATE=${WARN_MIN_ENRICHMENT_MAPPING_RATE:-0.30}
: "${PROJECT_DIR:?PROJECT_DIR is required}" "${SAMPLE_TABLE:?SAMPLE_TABLE is required}"
: "${MIN_COUNT:?MIN_COUNT is required}" "${MIN_SAMPLES:?MIN_SAMPLES is required}"
: "${FDR:?FDR is required}" "${LFC_CUTOFF:?LFC_CUTOFF is required}"
: "${CONTROL_GROUP:?CONTROL_GROUP is required}" "${TREATMENT_GROUP:?TREATMENT_GROUP is required}"
mkdir -p "$PROJECT_DIR/results" "$PROJECT_DIR/figures"

need() { command -v "$1" >/dev/null 2>&1 || { echo "Missing command: $1" >&2; exit 1; }; }
[[ -s "$SAMPLE_TABLE" ]] || { echo "Missing sample table: $SAMPLE_TABLE" >&2; exit 1; }

if [[ "$MODE" == full ]]; then
  : "${RAW_DIR:?RAW_DIR is required}" "${GENOME_FASTA:?GENOME_FASTA is required}" "${GTF_FILE:?GTF_FILE is required}"
  : "${HISAT2_INDEX_PREFIX:?HISAT2_INDEX_PREFIX is required}" "${THREADS:?THREADS is required}" "${STRANDED:?STRANDED is required}"
  need fastp; need hisat2; need hisat2-build; need samtools; need featureCounts
  [[ -s "$GENOME_FASTA" && -s "$GTF_FILE" ]] || { echo "Missing reference FASTA or GTF" >&2; exit 1; }
  if [[ "${BUILD_INDEX:-1}" == 1 && ! -s "${HISAT2_INDEX_PREFIX}.1.ht2" ]]; then
    mkdir -p "$(dirname "$HISAT2_INDEX_PREFIX")"
    ss="${HISAT2_INDEX_PREFIX}.ss"; exon="${HISAT2_INDEX_PREFIX}.exon"
    if command -v hisat2_extract_splice_sites.py >/dev/null 2>&1; then ss_tool=hisat2_extract_splice_sites.py; else ss_tool=extract_splice_sites.py; fi
    if command -v hisat2_extract_exons.py >/dev/null 2>&1; then exon_tool=hisat2_extract_exons.py; else exon_tool=extract_exons.py; fi
    need "$ss_tool"; need "$exon_tool"
    "$ss_tool" "$GTF_FILE" > "$ss"; "$exon_tool" "$GTF_FILE" > "$exon"
    hisat2-build -p "$THREADS" --ss "$ss" --exon "$exon" "$GENOME_FASTA" "$HISAT2_INDEX_PREFIX"
  fi
  [[ -s "${HISAT2_INDEX_PREFIX}.1.ht2" || -s "${HISAT2_INDEX_PREFIX}.1.ht2l" ]] || { echo "HISAT2 index is incomplete" >&2; exit 1; }
  clean="$PROJECT_DIR/results/clean_fastq"; bamdir="$PROJECT_DIR/results/alignment"
  mkdir -p "$clean/reports" "$bamdir"
  for sample in "${SAMPLES[@]}"; do
    r1="$RAW_DIR/${sample}${R1_SUFFIX}"; r2="$RAW_DIR/${sample}${R2_SUFFIX}"
    [[ -s "$r1" && -s "$r2" ]] || { echo "Missing FASTQ pair for $sample" >&2; exit 1; }
    fastp --in1 "$r1" --in2 "$r2" --out1 "$clean/${sample}_1.clean.fastq.gz" --out2 "$clean/${sample}_2.clean.fastq.gz" \
      --html "$clean/reports/${sample}_fastp.html" --json "$clean/reports/${sample}_fastp.json" --thread "$THREADS" \
      --qualified_quality_phred 20 --unqualified_percent_limit 40 --length_required 50 --cut_front --cut_tail \
      --cut_window_size 4 --cut_mean_quality 20 --detect_adapter_for_pe --correction
    hisat2 -x "$HISAT2_INDEX_PREFIX" -1 "$clean/${sample}_1.clean.fastq.gz" -2 "$clean/${sample}_2.clean.fastq.gz" \
      -p "$THREADS" --dta --summary-file "$bamdir/${sample}_alignment_summary.txt" | samtools sort -@ "$THREADS" -o "$bamdir/${sample}.sorted.bam" -
    samtools index -@ "$THREADS" "$bamdir/${sample}.sorted.bam"
  done
  bams=(); for sample in "${SAMPLES[@]}"; do bams+=("$bamdir/${sample}.sorted.bam"); done
  countdir="$PROJECT_DIR/results/counts"; mkdir -p "$countdir"
  featureCounts -a "$GTF_FILE" -o "$countdir/gene_counts.txt" -p --countReadPairs -T "$THREADS" -t exon -g gene_id -s "$STRANDED" "${bams[@]}"
  awk 'BEGIN{OFS="\t"} /^#/ {next} !seen {printf "Geneid"; for(i=7;i<=NF;i++){x=$i;sub(/^.*\//,"",x);sub(/\.sorted\.bam$/,"",x);printf OFS x};print "";seen=1;next} {printf $1;for(i=7;i<=NF;i++)printf OFS $i;print ""}' "$countdir/gene_counts.txt" > "$countdir/gene_counts_matrix.tsv"
  COUNTS_MATRIX="$countdir/gene_counts_matrix.tsv"
fi

need python3; need Rscript
[[ -s "${COUNTS_MATRIX:-}" ]] || { echo "Missing count matrix: ${COUNTS_MATRIX:-unset}" >&2; exit 1; }
filtered="$PROJECT_DIR/results/gene_counts_matrix_filtered.tsv"
python3 "$SKILL_DIR/filter_counts.py" --input "$COUNTS_MATRIX" --output "$filtered" --min-count "$MIN_COUNT" --min-samples "$MIN_SAMPLES"
Rscript "$SKILL_DIR/deseq2_report.R" --counts "$filtered" --samples "$SAMPLE_TABLE" --outdir "$PROJECT_DIR" --control "$CONTROL_GROUP" --treatment "$TREATMENT_GROUP" --fdr "$FDR" --lfc "$LFC_CUTOFF"

if [[ "${RUN_ENRICHMENT:-0}" == 1 ]]; then
  for f in "$ANNOTATION_TSV" "$GO_OBO" "$KO_TO_PATHWAY" "$PATHWAY_NAMES"; do [[ -s "$f" ]] || { echo "Missing enrichment input: $f" >&2; exit 1; }; done
  python3 "$SKILL_DIR/enrich_go_kegg.py" --annotation "$ANNOTATION_TSV" --go-obo "$GO_OBO" --ko-pathway "$KO_TO_PATHWAY" --pathway-names "$PATHWAY_NAMES" --background "$PROJECT_DIR/results/deseq2_results.tsv" --up "$PROJECT_DIR/results/deseq2_up_regulated_genes.tsv" --down "$PROJECT_DIR/results/deseq2_down_regulated_genes.tsv" --outdir "$PROJECT_DIR/results" --gene-col "$ANNOT_GENE_COL" --ko-col "$ANNOT_KO_COL" --go-col "$ANNOT_GO_COL"
  Rscript "$SKILL_DIR/plot_enrichment.R" --results "$PROJECT_DIR/results" --figures "$PROJECT_DIR/figures" --top-n "${TOP_ENRICHMENT_TERMS:-10}"
fi
manifest="$PROJECT_DIR/results/run_manifest.tsv"
{ printf 'key\tvalue\n'; printf 'run_date\t%s\n' "$(date '+%F %T %z')"; printf 'analysis_mode\t%s\n' "$MODE"; printf 'project_name\t%s\n' "$PROJECT_NAME"; printf 'species\t%s\n' "$SPECIES_NAME"; printf 'reference_fasta\t%s\n' "${GENOME_FASTA:-not used}"; printf 'annotation_gtf\t%s\n' "${GTF_FILE:-not used}"; printf 'reference_version\t%s\n' "$REFERENCE_VERSION"; printf 'config_file\t%s\n' "$CONFIG"; printf 'sample_table\t%s\n' "$SAMPLE_TABLE"; printf 'threads\t%s\n' "${THREADS:-not used}"; printf 'strandedness\t%s\n' "${STRANDED:-not used}"; printf 'min_count\t%s\n' "$MIN_COUNT"; printf 'min_samples\t%s\n' "$MIN_SAMPLES"; printf 'fdr\t%s\n' "$FDR"; printf 'lfc_cutoff\t%s\n' "$LFC_CUTOFF"; for tool_name in fastp hisat2 samtools featureCounts Rscript python3; do if command -v "$tool_name" >/dev/null 2>&1; then if [[ "$tool_name" == featureCounts ]]; then version=$(featureCounts -h 2>&1 | awk '/^Version / {print; exit}' || true); else version=$($tool_name --version 2>&1 | head -1 || true); fi; printf 'software_%s\t%s\n' "$tool_name" "${version:-available}"; fi; done; } > "$manifest"
python3 "$SKILL_DIR/generate_rnaseq_report.py" --project "$PROJECT_DIR" --samples "$SAMPLE_TABLE" --control "$CONTROL_GROUP" --treatment "$TREATMENT_GROUP" --output "$PROJECT_DIR/RNA-seq分析报告.md" --project-name "$PROJECT_NAME" --species "$SPECIES_NAME" --mode "$MODE" --config "$CONFIG" --reference "${GENOME_FASTA:-not used}" --annotation "${GTF_FILE:-not used}" --reference-version "$REFERENCE_VERSION" --strandedness "${STRANDED:-not used}" --min-count "$MIN_COUNT" --min-samples "$MIN_SAMPLES" --fdr "$FDR" --lfc "$LFC_CUTOFF" --warn-min-replicates "$WARN_MIN_REPLICATES" --warn-min-retention "$WARN_MIN_READ_RETENTION" --warn-min-q30 "$WARN_MIN_Q30" --warn-min-alignment "$WARN_MIN_ALIGNMENT_RATE" --warn-min-assignment "$WARN_MIN_ASSIGNMENT_RATE" --warn-min-filtered-fraction "$WARN_MIN_FILTERED_GENE_FRACTION" --warn-max-depth-ratio "$WARN_MAX_GROUP_LIBRARY_DEPTH_RATIO" --warn-min-within-correlation "$WARN_MIN_WITHIN_GROUP_CORRELATION" --warn-min-mapping-rate "$WARN_MIN_ENRICHMENT_MAPPING_RATE"
echo "Completed. Results: $PROJECT_DIR/results ; figures: $PROJECT_DIR/figures"
