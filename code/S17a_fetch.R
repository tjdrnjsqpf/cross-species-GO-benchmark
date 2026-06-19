#!/usr/bin/env Rscript
# S17a fetch — for each Bgee target species, download RNA-seq via BgeeDB and aggregate to a slim
# gene x anatomy mean-rank table (Bgee Rank = cross-condition normalised expression; lower=higher).
# Writes data/bgee/<bgee_name>.slim.tsv. Skips species already done.
suppressMessages({library(BgeeDB); library(data.table)})
ROOT <- "/var2/lsg/Claude_Code/Cross-species-GeneOntology"
setwd(file.path(ROOT,"data/bgee"))
tg <- fread(file.path(ROOT,"data/bgee/targets.tsv"))
species <- unique(tg$bgee_name)
for (sp in species) {
  out <- paste0(sp, ".slim.tsv")
  if (file.exists(out) && file.info(out)$size > 1000) { cat("[done]", sp, "\n"); next }
  cat("[fetch]", sp, "\n")
  ok <- tryCatch({
    bg <- Bgee$new(species=sp, dataType="rna_seq")
    d <- suppressWarnings(getData(bg))
    if (is.list(d) && !is.data.frame(d)) d <- rbindlist(d, fill=TRUE)
    setDT(d)
    d <- d[!is.na(Rank) & Rank>0, .(rank=mean(Rank), tpm=mean(TPM, na.rm=TRUE), n=.N),
           by=.(Gene.ID, Anatomical.entity.ID)]
    fwrite(d, out, sep="\t")
    cat("  ok", sp, "rows=", nrow(d), "\n"); TRUE
  }, error=function(e){cat("  FAIL", sp, conditionMessage(e), "\n"); FALSE})
}
cat("S17A FETCH DONE\n")
