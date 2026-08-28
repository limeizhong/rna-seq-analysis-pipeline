---
name: rna-seq-analysis-pipeline
description: Analyze conventional paired-end bulk RNA-seq data reproducibly, from FASTQ quality control, HISAT2 genome alignment and featureCounts gene quantification through DESeq2 differential expression, QC plots, DEG plots, and optional local GO/KEGG enrichment. Use when asked to process RNA-seq FASTQ files, create a gene-count matrix, compare treatment vs control RNA-seq groups, run DESeq2, or migrate/re-run this workflow in another project or GitHub clone.
---

# RNA-seq analysis pipeline

Use this skill for conventional, paired-end, reference-genome RNA-seq. It reproduces the validated workflow from experiment08 and experiment09: fastp → HISAT2 → featureCounts → DESeq2/QC → optional GO/KEGG enrichment.

## Start a project

Copy the bundled configuration and edit paths, sample names, groups, reference files, and thresholds. Keep raw FASTQ and reference files outside the skill directory; put only configuration and generated results in the project directory.

```bash
cp <skill-dir>/scripts/rnaseq_config.example.sh ./rnaseq_config.sh
# edit rnaseq_config.sh
bash <skill-dir>/scripts/run_rnaseq.sh --config ./rnaseq_config.sh --mode full
```

Use `--mode counts` when a raw integer gene-count matrix already exists. The count matrix must be tab-delimited, genes in the first column, samples in remaining columns. `sample_group.tsv` must contain `sample` and `group` columns; its sample names must exactly match the count-matrix columns.

## Required checks

- Confirm that the input is **raw integer gene counts**, not TPM/FPKM/log-transformed values, before DESeq2.
- Set `STRANDED` from the library preparation; never assume unstranded (`0`) without checking.
- Set `CONTROL_GROUP` and `TREATMENT_GROUP` explicitly. Positive log2FoldChange always means treatment/control.
- Keep at least two biological replicates per group; three or more are strongly preferred.
- Review fastp reports, alignment summaries, assignment rates, PCA, correlation heatmap, and clustering before interpreting DEGs. Flag outliers rather than silently removing them.

## Report and quality review

Treat `RNA-seq分析报告.md` as a scientific result record, not a decorative summary. The runner collects structured QC, alignment, count-assignment, VST/PCA/correlation/clustering, DESeq2, enrichment-mapping, software-version, and configuration records before rendering the report. It separates observed values, empirical quality warnings, and cautious biological interpretation.

Set the `WARN_*` values in the project configuration to make empirical warnings appropriate for the library type and study. Do not interpret a warning as proof of a cause, or delete a sample solely from PCA or clustering. Read the warning together with raw QC reports and experimental metadata.

## Outputs

The runner writes a stable `results/` and `figures/` layout: cleaned FASTQ, BAM/BAI and HISAT2 summaries, `gene_counts_matrix.tsv`, filtered counts, DESeq2 tables, DEG tables, VST matrices, and PCA/correlation/cluster/MA/volcano/DEG-heatmap figures. Optional enrichment additionally writes separate up/down GO and KEGG tables plus bar plots and bubble plots for the top terms in each direction. It always writes `RNA-seq分析报告.md` at the project root, covering quality control, alignment, quantification, sample QC, differential expression, and functional enrichment.

Use the filtered genes that actually enter DESeq2 as the enrichment background, not the entire genome. Treat cross-species KEGG pathway labels cautiously.

## Resources

- `scripts/rnaseq_config.example.sh`: portable project configuration template.
- `scripts/run_rnaseq.sh`: one-command full or count-matrix workflow.
- `scripts/filter_counts.py`: low-expression filter.
- `scripts/deseq2_report.R`: DESeq2, sample QC, classification, and figures.
- `scripts/enrich_go_kegg.py`: optional local GO/KEGG hypergeometric enrichment.
- `scripts/plot_enrichment.R`: enrichment bar plots and bubble plots.
- `scripts/generate_rnaseq_report.py`: final Chinese Markdown report with quality warnings, figure captions, interpretation, limitations, and reproducibility record.
- `references/input-and-environment.md`: input layout, dependencies, and implementation notes.
