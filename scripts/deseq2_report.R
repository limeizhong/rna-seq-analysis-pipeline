#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(DESeq2); library(ggplot2); library(pheatmap) })

opts <- commandArgs(trailingOnly=TRUE)
arg <- function(name) { i <- match(name, opts); if (is.na(i) || i == length(opts)) stop("Missing ", name); opts[i + 1] }
count_file <- arg("--counts"); sample_file <- arg("--samples"); outdir <- arg("--outdir")
control <- arg("--control"); treatment <- arg("--treatment")
fdr <- as.numeric(arg("--fdr")); lfc_cut <- as.numeric(arg("--lfc"))
resdir <- file.path(outdir, "results"); figdir <- file.path(outdir, "figures")
dir.create(resdir, recursive=TRUE, showWarnings=FALSE); dir.create(figdir, recursive=TRUE, showWarnings=FALSE)

counts <- as.matrix(read.table(count_file, header=TRUE, row.names=1, sep="\t", check.names=FALSE, quote=""))
if (anyNA(counts) || any(counts < 0) || any(counts != round(counts))) stop("Counts must be non-negative integers")
storage.mode(counts) <- "integer"
sample_info <- read.table(sample_file, header=TRUE, sep="", stringsAsFactors=FALSE, check.names=FALSE, quote="", comment.char="", fill=TRUE)
if (!all(c("sample", "group") %in% names(sample_info))) stop("Sample table must have sample and group columns")
sample_info$sample <- trimws(sample_info$sample); sample_info$group <- trimws(sample_info$group)
if (anyDuplicated(sample_info$sample) || !setequal(colnames(counts), sample_info$sample)) stop("Count-matrix sample names must exactly match sample table")
sample_info <- sample_info[match(colnames(counts), sample_info$sample),,drop=FALSE]
if (!all(c(control, treatment) %in% sample_info$group)) stop("Control/treatment group is absent from sample table")
sample_info$group <- factor(sample_info$group, levels=c(control, treatment)); rownames(sample_info) <- sample_info$sample

dds <- DESeqDataSetFromMatrix(countData=counts, colData=sample_info, design=~group)
dds <- DESeq(dds)
res <- results(dds, contrast=c("group", treatment, control), alpha=fdr)
res_df <- data.frame(Geneid=rownames(res), as.data.frame(res), check.names=FALSE)
res_df <- res_df[order(res_df$pvalue, na.last=TRUE),]
reg <- rep("Not significant", nrow(res_df)); ok <- !is.na(res_df$padj) & res_df$padj < fdr
reg[ok & res_df$log2FoldChange >= lfc_cut] <- "Up-regulated"
reg[ok & res_df$log2FoldChange <= -lfc_cut] <- "Down-regulated"
res_df$regulation <- reg
write.table(res_df[, names(res_df) != "regulation"], file.path(resdir,"deseq2_results.tsv"), sep="\t", quote=FALSE, row.names=FALSE, na="NA")
write.table(res_df, file.path(resdir,"deseq2_genes_classified.tsv"), sep="\t", quote=FALSE, row.names=FALSE, na="NA")
write.table(res_df[reg != "Not significant",], file.path(resdir,"deseq2_significant_genes.tsv"), sep="\t", quote=FALSE, row.names=FALSE, na="NA")
write.table(res_df[reg == "Up-regulated",], file.path(resdir,"deseq2_up_regulated_genes.tsv"), sep="\t", quote=FALSE, row.names=FALSE, na="NA")
write.table(res_df[reg == "Down-regulated",], file.path(resdir,"deseq2_down_regulated_genes.tsv"), sep="\t", quote=FALSE, row.names=FALSE, na="NA")
norm <- as.data.frame(counts(dds, normalized=TRUE)); norm$Geneid <- rownames(norm); norm <- norm[,c("Geneid", colnames(counts))]
write.table(norm, file.path(resdir,"deseq2_normalized_counts.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# Blind VST is used only for visualization; DESeq2 testing remains on raw counts.
vsd <- varianceStabilizingTransformation(dds, blind=TRUE); vst <- assay(vsd)
write.table(data.frame(Geneid=rownames(vst), vst, check.names=FALSE), file.path(resdir,"vst_expression_matrix.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
pca <- prcomp(t(vst)); pct <- round(100*pca$sdev^2/sum(pca$sdev^2), 1)
pca_df <- data.frame(sample=rownames(pca$x), PC1=pca$x[,1], PC2=pca$x[,2], group=sample_info[rownames(pca$x),"group"])
pca_df$PC1_variance_percent <- pct[1]; pca_df$PC2_variance_percent <- pct[2]
write.table(pca_df, file.path(resdir,"pca_coordinates.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cols <- c(control="#2a78d6", treatment="#e34948"); names(cols) <- c(control,treatment)
p <- ggplot(pca_df,aes(PC1,PC2,color=group,label=sample))+geom_point(size=3)+geom_text(vjust=-0.7)+scale_color_manual(values=cols)+theme_bw()+labs(title="PCA (VST)",x=paste0("PC1: ",pct[1],"%"),y=paste0("PC2: ",pct[2],"%"))
ggsave(file.path(figdir,"PCA_plot.png"),p,width=7,height=6,dpi=300); ggsave(file.path(figdir,"PCA_plot.pdf"),p,width=7,height=6)
cor_mat <- cor(vst, method="pearson")
write.table(data.frame(sample=rownames(cor_mat), cor_mat, check.names=FALSE), file.path(resdir,"sample_correlation_matrix.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
pheatmap(cor_mat, filename=file.path(figdir,"correlation_heatmap.png"), width=7, height=6, main="Sample Pearson correlation")
pheatmap(cor_mat, filename=file.path(figdir,"correlation_heatmap.pdf"), width=7, height=6, main="Sample Pearson correlation")
sample_hc <- hclust(dist(t(vst)), method="complete")
clusters <- cutree(sample_hc, k=min(2,ncol(vst)))
sample_qc <- data.frame(sample=colnames(vst), group=sample_info[colnames(vst),"group"], cluster_k2=clusters[colnames(vst)], stringsAsFactors=FALSE)
sample_qc$mean_within_group_r <- sapply(seq_len(nrow(sample_qc)), function(i) { s <- sample_qc$sample[i]; x <- sample_qc$sample[sample_qc$group == sample_qc$group[i] & sample_qc$sample != s]; if(length(x)) mean(cor_mat[s,x]) else NA })
sample_qc$mean_between_group_r <- sapply(seq_len(nrow(sample_qc)), function(i) { s <- sample_qc$sample[i]; x <- sample_qc$sample[sample_qc$group != sample_qc$group[i]]; if(length(x)) mean(cor_mat[s,x]) else NA })
write.table(sample_qc, file.path(resdir,"sample_qc_metrics.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
png(file.path(figdir,"hierarchical_clustering.png"), width=2100,height=1600,res=300); plot(sample_hc, main="Sample clustering (VST, Euclidean distance, complete linkage)"); dev.off()
pdf(file.path(figdir,"hierarchical_clustering.pdf"), width=7,height=5.5); plot(sample_hc, main="Sample clustering (VST, Euclidean distance, complete linkage)"); dev.off()

plot_df <- res_df; plot_df$negLog10Padj <- -log10(pmax(ifelse(is.na(plot_df$padj),1,plot_df$padj), .Machine$double.xmin)); plot_df$log10BaseMean <- log10(plot_df$baseMean+1)
write.table(plot_df, file.path(resdir,"deseq2_plot_data.tsv"), sep="\t", quote=FALSE,row.names=FALSE,na="NA")
colors <- c("Up-regulated"="#e34948","Down-regulated"="#2a78d6","Not significant"="#b6b6b6")
ma <- ggplot(plot_df,aes(log10BaseMean,log2FoldChange,color=regulation))+geom_point(alpha=.7)+geom_hline(yintercept=c(-lfc_cut,lfc_cut),linetype="dashed")+scale_color_manual(values=colors)+theme_bw()+labs(title="MA plot",x="log10(baseMean + 1)")
vol <- ggplot(plot_df,aes(log2FoldChange,negLog10Padj,color=regulation))+geom_point(alpha=.7)+geom_hline(yintercept=-log10(fdr),linetype="dashed")+geom_vline(xintercept=c(-lfc_cut,lfc_cut),linetype="dashed")+scale_color_manual(values=colors)+theme_bw()+labs(title="Volcano plot",y="-log10(padj)")
for (x in list(list(name="MA_plot", plot=ma), list(name="volcano_plot", plot=vol))) { ggsave(file.path(figdir,paste0(x$name,".png")),x$plot,width=7,height=6,dpi=300); ggsave(file.path(figdir,paste0(x$name,".pdf")),x$plot,width=7,height=6) }

deg <- res_df$Geneid[reg != "Not significant"]
if (length(deg) > 1) {
  z <- t(scale(t(vst[deg,,drop=FALSE]))); z[is.na(z)] <- 0
  write.table(data.frame(Geneid=rownames(z), z, check.names=FALSE), file.path(resdir,"deseq2_deg_zscore_matrix.tsv"), sep="\t",quote=FALSE,row.names=FALSE)
  ann <- data.frame(Group=sample_info$group,row.names=sample_info$sample)
  pheatmap(z, annotation_col=ann, filename=file.path(figdir,"heatmap_deg.png"),width=7,height=max(5,min(12,nrow(z)*.18)),main="DEG heatmap (VST row z-score)")
  pheatmap(z, annotation_col=ann, filename=file.path(figdir,"heatmap_deg.pdf"),width=7,height=max(5,min(12,nrow(z)*.18)),main="DEG heatmap (VST row z-score)")
}
summary <- c(sprintf("DESeq2 version: %s", as.character(packageVersion("DESeq2"))),"Design formula: ~ group",sprintf("Comparison: %s vs %s",treatment,control),"Test: Wald test; adjustment: Benjamini-Hochberg; LFC shrinkage: not applied",sprintf("Genes entering DESeq2: %d",nrow(dds)),sprintf("Genes with non-missing pvalue: %d",sum(!is.na(res_df$pvalue))),sprintf("Genes with non-missing padj: %d",sum(!is.na(res_df$padj))),sprintf("Independent-filtered/outlier pvalue NA: %d",sum(is.na(res_df$pvalue))),sprintf("Significant (padj < %g): %d",fdr,sum(ok)),sprintf("DEG (padj < %g and |log2FC| >= %g): %d",fdr,lfc_cut,length(deg)),sprintf("Up: %d; Down: %d",sum(reg=="Up-regulated"),sum(reg=="Down-regulated")))
writeLines(summary,file.path(resdir,"deseq2_summary.txt")); cat(paste(summary,collapse="\n"),"\n")
