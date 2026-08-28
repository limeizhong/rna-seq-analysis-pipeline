# Inputs and environment

Use a Conda environment containing `fastp`, `hisat2`, `samtools`, `subread` (featureCounts), `R`, and `python3`. Install R/Bioconductor packages `DESeq2`, `ggplot2`, and `pheatmap`.

The full mode expects `RAW_DIR/<sample>_1.fastq.gz` and `<sample>_2.fastq.gz`; change the suffix settings if needed. The GTF must include exon records and a valid `gene_id` attribute. `featureCounts --countReadPairs` is intentional: for paired-end RNA-seq, it prevents version-dependent read-level rather than fragment-level counts.

`sample_group.tsv` example (tab or whitespace delimited):

```text
sample	group
sample1	Treatment
sample2	Treatment
sample3	Treatment
sample4	Control
sample5	Control
sample6	Control
```

`--mode counts` is appropriate when quantification is already complete. It still produces low-expression filtering, DESeq2, PCA/correlation/clustering, MA/volcano, and DEG heatmap outputs.

Optional enrichment is intentionally local and annotation-version-dependent. Its default annotation-column indices match the soybean `info_annot.txt` used in experiment09. Change the three indices to match another annotation file, and check KO-to-pathway cross-species labels before biological interpretation.

When `RUN_ENRICHMENT=1`, the runner also produces `GO_enrichment_barplot`, `KEGG_enrichment_barplot`, `GO_enrichment_bubbleplot`, and `KEGG_enrichment_bubbleplot` (PNG and PDF) in `figures/`. Terms are selected by nominal p-value for visibility; use the asterisk and reported `padj`, not rank alone, to determine FDR significance.

Every completed run writes `RNA-seq分析报告.md` to the project root. It contains these fixed sections: sequencing-data quality control, sequence alignment, gene-expression quantification, sample quality evaluation, differential expression analysis, and functional enrichment analysis. It links or embeds every figure that exists; when a count-matrix-only run lacks FASTQ-derived metrics, the report explicitly marks those values unavailable.

Configure `PROJECT_NAME`, `SPECIES_NAME`, and `REFERENCE_VERSION` for a complete project record. The `WARN_*` settings govern empirical report alerts for replicate count, read retention, Q30, alignment and assignment rate, retained-gene fraction, library-depth imbalance, within-group correlation, and annotation mapping rate. They are adjustable warnings, not universal acceptance criteria.
