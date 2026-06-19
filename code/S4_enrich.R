#!/usr/bin/env Rscript
# S4 — enrichment for every (annotation set x gene list x aspect) via clusterProfiler::enricher.
# Annotation sets: truth (focal EXP) + all transferred[ref,method,evset].
# Background = focal measured proteome (same for all -> fair set comparison; Primmer Box 3).
suppressMessages({library(clusterProfiler); library(data.table)})

ROOT <- "/var2/lsg/Claude_Code/Cross-species-GeneOntology"
args <- commandArgs(trailingOnly=TRUE)
# args: manifest(label\tpath) lists.tsv background.tsv out.tsv
manifest <- fread(args[1], header=TRUE)
lists    <- fread(args[2], header=TRUE)        # list_id, focal_acc
bg       <- fread(args[3], header=TRUE)[[1]]
outpath  <- args[4]

NS <- c(biological_process="BP", molecular_function="MF", cellular_component="CC")
list_ids <- unique(lists$list_id)
res_all <- list(); k <- 0

for (i in seq_len(nrow(manifest))) {
  label <- manifest$label[i]; path <- manifest$path[i]
  if (!file.exists(path)) { message("missing ", path); next }
  ann <- fread(path, header=TRUE)             # col1 gene, col2 go, col3 namespace
  setnames(ann, 1:3, c("gene","go","ns"))
  for (nsname in names(NS)) {
    asp <- NS[[nsname]]
    sub <- ann[ns==nsname]
    if (nrow(sub)==0) next
    t2g <- sub[, .(go, gene)]
    for (lid in list_ids) {
      genes <- lists[list_id==lid, focal_acc]
      genes <- intersect(genes, bg)
      if (length(genes) < 5) next
      e <- tryCatch(enricher(gene=genes, universe=bg, TERM2GENE=t2g,
                pvalueCutoff=1, qvalueCutoff=1, minGSSize=5, maxGSSize=2000,
                pAdjustMethod="BH"), error=function(x) NULL)
      if (is.null(e) || nrow(as.data.frame(e))==0) next
      df <- as.data.table(as.data.frame(e))
      df <- df[pvalue < 0.1]                    # bound output; keeps sig + near-sig
      if (nrow(df)==0) next
      df[, rank := frank(pvalue, ties.method="average")]
      k <- k+1
      res_all[[k]] <- data.table(label=label, list_id=lid, aspect=asp,
                                 go_id=df$ID, pvalue=df$pvalue, padj=df$p.adjust, rank=df$rank)
    }
  }
  message(sprintf("[S4] done %s", label))
}

out <- rbindlist(res_all)
fwrite(out, outpath, sep="\t")
cat(sprintf("[S4] wrote %d rows -> %s\n", nrow(out), outpath))
