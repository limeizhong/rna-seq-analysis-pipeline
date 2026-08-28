#!/usr/bin/env bash
# Copy this file into an RNA-seq project and edit every value marked CHANGE ME.
# Paths may be absolute or relative to the directory from which run_rnaseq.sh is called.

# Required project inputs
PROJECT_DIR="$(pwd)"
PROJECT_NAME="RNA-seq project"                 # shown in final report
SPECIES_NAME="Unknown species"                  # e.g. Glycine max
REFERENCE_VERSION="Not specified"               # assembly/annotation release or source
RAW_DIR="${PROJECT_DIR}/data/fastq"                 # CHANGE ME; ignored with --mode counts
COUNTS_MATRIX="${PROJECT_DIR}/data/gene_counts_matrix.tsv"  # CHANGE ME; ignored with --mode full
SAMPLE_TABLE="${PROJECT_DIR}/data/sample_group.tsv" # columns: sample<TAB>group
SAMPLES=(sample1 sample2 sample3 sample4 sample5 sample6)    # CHANGE ME

# Paired-end FASTQ naming: RAW_DIR/<sample>${R1_SUFFIX} and ...R2_SUFFIX
R1_SUFFIX="_1.fastq.gz"
R2_SUFFIX="_2.fastq.gz"

# Required only with --mode full
GENOME_FASTA="${PROJECT_DIR}/reference/genome.fa"  # CHANGE ME
GTF_FILE="${PROJECT_DIR}/reference/annotation.gtf" # exon records must contain gene_id
HISAT2_INDEX_PREFIX="${PROJECT_DIR}/reference/hisat2_index/genome"
BUILD_INDEX=1                                       # set 0 if index files already exist

# Analysis settings
THREADS=8
STRANDED=0              # featureCounts: 0 unstranded, 1 stranded, 2 reversely stranded
MIN_COUNT=1             # teaching-data default; raise for deeper libraries if justified
MIN_SAMPLES=3           # normally the smallest biological-group size
FDR=0.05
LFC_CUTOFF=1.0
CONTROL_GROUP="Control"     # must exactly match sample_group.tsv
TREATMENT_GROUP="Treatment" # positive log2FC = treatment/control

# Report-only empirical warning thresholds; adjust for the project and interpret as alerts, not hard rules.
WARN_MIN_REPLICATES=3
WARN_MIN_READ_RETENTION=70
WARN_MIN_Q30=80
WARN_MIN_ALIGNMENT_RATE=70
WARN_MIN_ASSIGNMENT_RATE=50
WARN_MIN_FILTERED_GENE_FRACTION=0.05
WARN_MAX_GROUP_LIBRARY_DEPTH_RATIO=2.0
WARN_MIN_WITHIN_GROUP_CORRELATION=0.80
WARN_MIN_ENRICHMENT_MAPPING_RATE=0.30

# Optional local GO/KEGG enrichment. Set RUN_ENRICHMENT=1 only after supplying all paths.
RUN_ENRICHMENT=0
TOP_ENRICHMENT_TERMS=10   # most nominally enriched terms shown per direction in each plot
ANNOTATION_TSV="${PROJECT_DIR}/reference/annotation.tsv"  # original teaching file layout
GO_OBO="${PROJECT_DIR}/reference/go-basic.obo"
KO_TO_PATHWAY="${PROJECT_DIR}/reference/ko_to_pathway.tsv"
PATHWAY_NAMES="${PROJECT_DIR}/reference/kegg_pathway_names.tsv"
# Zero-based columns in ANNOTATION_TSV; the teaching soybean annotation uses gene=1, KO=8, GO=9.
ANNOT_GENE_COL=1
ANNOT_KO_COL=8
ANNOT_GO_COL=9
