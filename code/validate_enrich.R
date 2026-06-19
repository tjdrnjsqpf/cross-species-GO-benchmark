#!/usr/bin/env Rscript
# Validate fast hypergeometric enricher against clusterProfiler::enricher on truth annotation.
suppressMessages({library(clusterProfiler); library(data.table)})
ROOT <- Sys.getenv("GOTX_ROOT", unset = "..")
ann <- fread(file.path(ROOT,"results/truth/zebrafish_truth_annotation.tsv"))
setnames(ann,1:3,c("gene","go","ns"))
lists <- fread(file.path(ROOT,"results/enrichment/lists_synthetic.tsv"))
bg <- fread(file.path(ROOT,"results/enrichment/background.tsv"))[[1]]
sel <- head(unique(lists$list_id), 5)
NS <- c(biological_process="BP")  # validate BP
res <- list(); k<-0
for (lid in sel) {
  genes <- intersect(lists[list_id==lid, focal_acc], bg)
  sub <- ann[ns=="biological_process"]
  e <- enricher(genes, universe=bg, TERM2GENE=sub[,.(go,gene)],
                pvalueCutoff=1, qvalueCutoff=1, minGSSize=5, maxGSSize=2000, pAdjustMethod="BH")
  df <- as.data.table(as.data.frame(e))
  if (nrow(df)==0) next
  k<-k+1; res[[k]] <- data.table(list_id=lid, go_id=df$ID, cp_pvalue=df$pvalue, cp_padj=df$p.adjust)
}
fwrite(rbindlist(res), file.path(ROOT,"results/enrichment/validate_cp.tsv"), sep="\t")
cat("[validate] clusterProfiler ref written for", k, "lists\n")
