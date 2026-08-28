#!/usr/bin/env Rscript
# Draw bar and bubble plots from enrichment_GO_{up,down}.tsv and enrichment_KEGG_{up,down}.tsv.
suppressPackageStartupMessages(library(ggplot2))

opts <- commandArgs(trailingOnly=TRUE)
arg <- function(name, default=NULL) { i <- match(name,opts); if (is.na(i)) return(default); if (i==length(opts)) stop("Missing value for ",name); opts[i+1] }
resdir <- arg("--results"); figdir <- arg("--figures"); top_n <- as.integer(arg("--top-n", "10"))
dir.create(figdir, recursive=TRUE, showWarnings=FALSE)

read_enrichment <- function(kind, direction) {
  file <- file.path(resdir, sprintf("enrichment_%s_%s.tsv", kind, tolower(direction)))
  if (!file.exists(file) || file.info(file)$size == 0) return(NULL)
  d <- read.table(file, header=TRUE, sep="\t", quote="", check.names=FALSE, stringsAsFactors=FALSE)
  if (!nrow(d)) return(NULL)
  d$Direction <- direction
  d$TermLabel <- if (kind == "GO") paste0(d$Term, " (", d$ontology, ")") else d$Pathway
  d$TermLabel <- vapply(d$TermLabel, function(x) paste(strwrap(x, width=46), collapse="\n"), character(1))
  ratio <- strsplit(d$GeneRatio, "/", fixed=TRUE)
  d$GeneRatioValue <- vapply(ratio, function(x) as.numeric(x[1])/as.numeric(x[2]), numeric(1))
  d$Count <- vapply(ratio, function(x) as.numeric(x[1]), numeric(1))
  d$negLog10P <- -log10(pmax(d$pvalue, .Machine$double.xmin))
  d$negLog10Padj <- -log10(pmax(d$padj, .Machine$double.xmin))
  d$FDR <- d$padj < 0.05
  d
}

top_terms <- function(d) {
  if (is.null(d)) return(NULL)
  by_direction <- split(d, d$Direction)
  do.call(rbind, lapply(by_direction, function(x) x[order(x$pvalue),,drop=FALSE][seq_len(min(top_n,nrow(x))),,drop=FALSE]))
}

save_plot <- function(p, stem, height) {
  ggsave(file.path(figdir,paste0(stem,".png")),p,width=10,height=height,dpi=300)
  ggsave(file.path(figdir,paste0(stem,".pdf")),p,width=10,height=height)
}

make_plots <- function(kind) {
  d <- top_terms(do.call(rbind, Filter(Negate(is.null), list(read_enrichment(kind,"Up"),read_enrichment(kind,"Down")))))
  if (is.null(d) || !nrow(d)) { message("No ",kind," enrichment rows; no plot produced."); return(invisible()) }
  # Add direction to factor labels so the two facets can be independently ordered.
  d$facet_label <- paste(d$TermLabel,d$Direction,sep="___")
  d <- d[order(d$Direction,d$pvalue,decreasing=TRUE),,drop=FALSE]
  d$facet_label <- factor(d$facet_label, levels=unique(d$facet_label))
  strip_direction <- function(x) sub("___(Up|Down)$", "", x)
  height <- max(5, 2.3 + 0.42*nrow(d))
  theme_report <- theme_bw(base_size=12) + theme(panel.grid.minor=element_blank(), panel.grid.major.y=element_blank(), strip.background=element_rect(fill="grey92"), strip.text=element_text(face="bold"), axis.text.y=element_text(size=9), plot.title=element_text(face="bold",hjust=.5))
  if (kind == "GO") {
    fill <- scale_fill_manual(values=c(BP="#0072B2",MF="#D55E00",CC="#009E73"),name="Ontology")
    bar <- ggplot(d,aes(negLog10P,facet_label))+geom_col(aes(fill=ontology),width=.72)+geom_vline(xintercept=-log10(.05),linetype="dashed")+geom_text(data=d[d$FDR,,drop=FALSE],aes(label="*"),hjust=-.35,size=5)+fill
    bubble <- ggplot(d,aes(GeneRatioValue,facet_label))+geom_point(aes(size=Count,color=negLog10Padj),alpha=.9)+scale_color_viridis_c(name="-log10(padj)")+scale_size_continuous(name="DEG count")
  } else {
    fill <- scale_fill_manual(values=c(Up="#e34948",Down="#2a78d6"),name="Direction")
    bar <- ggplot(d,aes(negLog10P,facet_label))+geom_col(aes(fill=Direction),width=.72)+geom_vline(xintercept=-log10(.05),linetype="dashed")+geom_text(data=d[d$FDR,,drop=FALSE],aes(label="*"),hjust=-.35,size=5)+fill
    bubble <- ggplot(d,aes(GeneRatioValue,facet_label))+geom_point(aes(size=Count,color=negLog10Padj),alpha=.9)+scale_color_viridis_c(name="-log10(padj)")+scale_size_continuous(name="DEG count")
  }
  bar <- bar + facet_grid(Direction~.,scales="free_y",space="free_y") + scale_y_discrete(labels=strip_direction) + scale_x_continuous(expand=expansion(mult=c(0,.12))) + labs(title=paste(kind,"enrichment: top terms"),subtitle="Dashed line: nominal p = 0.05; *: padj < 0.05",x="-log10(pvalue)",y=NULL) + theme_report
  bubble <- bubble + facet_grid(Direction~.,scales="free_y",space="free_y") + scale_y_discrete(labels=strip_direction) + labs(title=paste(kind,"enrichment bubble plot: top terms"),x="Gene ratio",y=NULL) + theme_report
  save_plot(bar,paste0(kind,"_enrichment_barplot"),height)
  save_plot(bubble,paste0(kind,"_enrichment_bubbleplot"),height)
}

make_plots("GO")
make_plots("KEGG")
