#!/usr/bin/env Rscript
# S6 — deliverable figures + guide table from set-level metrics.
suppressMessages({library(data.table); library(ggplot2)})
ROOT <- Sys.getenv("GOTX_ROOT", unset = "..")
OUT <- Sys.getenv("GOTX_OUT", file.path(ROOT,"results"))
TRACK <- Sys.getenv("GOTX_TRACK", "fish")
FOCAL <- Sys.getenv("GOTX_FOCAL", "zebrafish")
EQD <- strsplit(Sys.getenv("GOTX_EQD", "chicken,mouse,human"), ",")[[1]]
TT <- function(x) paste0(x, "  [track: ", TRACK, ", focal: ", FOCAL, "]")
m <- fread(file.path(OUT,"metrics/setlevel_metrics.tsv"))
figdir <- file.path(OUT,"figures"); tabdir <- file.path(OUT,"tables")
dir.create(figdir, showWarnings=FALSE); dir.create(tabdir, showWarnings=FALSE)

asp_lab <- c(BP="BP (process)", MF="MF (function)", CC="CC (component)")
m[, aspect := factor(aspect, levels=c("MF","BP","CC"))]

## ---- Fig 1: divergence tolerance curve (all-IC), F1 vs %identity ----
d1 <- m[ic_bin=="all" & list_source=="synthetic", .(f1=mean(f1,na.rm=TRUE), recall=mean(recall,na.rm=TRUE),
                          se=sd(f1,na.rm=TRUE)/sqrt(.N)),
         by=.(ref,method,evset,aspect,median_pident,My)]
p1 <- ggplot(d1[evset=="all"], aes(median_pident, f1, color=aspect, group=aspect)) +
  geom_point(size=2) + geom_line() +
  geom_errorbar(aes(ymin=f1-se, ymax=f1+se), width=1, alpha=.5) +
  facet_wrap(~method, ncol=2) +
  scale_x_reverse() +  # closer (high identity) on left
  labs(title=TT("Divergence tolerance curve (all-evidence transfer)"),
       subtitle="x reversed: left=closer (high %identity). Mean set-level F1 +/- SE over synthetic lists",
       x="median ortholog % identity (focal-reference)", y="set-level F1 (transferred vs truth)") +
  theme_bw(base_size=11)
ggsave(file.path(figdir,"fig1_tolerance_curve.pdf"), p1, width=9, height=5)
ggsave(file.path(figdir,"fig1_tolerance_curve.svg"), p1, width=9, height=5)

## ---- Fig 1b: also vs My ----
p1b <- ggplot(d1[evset=="all"], aes(My, f1, color=aspect, group=aspect)) +
  geom_point(size=2) + geom_line() + facet_wrap(~method, ncol=2) +
  labs(title=TT("Divergence tolerance vs divergence time (all-evidence)"),
       x="divergence time from focal (My)", y="set-level F1") + theme_bw(base_size=11)
ggsave(file.path(figdir,"fig1b_tolerance_vs_My.pdf"), p1b, width=9, height=5)
ggsave(file.path(figdir,"fig1b_tolerance_vs_My.svg"), p1b, width=9, height=5)

## ---- Fig 2: distance vs richness decomposition (equidistant group) ----
eqd <- EQD
d2 <- m[ic_bin=="all" & list_source=="synthetic" & ref %in% eqd & method=="rbh",
        .(f1=mean(f1,na.rm=TRUE), recall=mean(recall,na.rm=TRUE), precision=mean(precision,na.rm=TRUE),
          richness=mean(ref_richness)), by=.(ref,evset,aspect)]
d2[, ref := factor(ref, levels=EQD)]  # low->high richness, ~equidistant
p2 <- ggplot(d2, aes(ref, recall, fill=evset)) +
  geom_col(position=position_dodge2(preserve="single"), width=0.8) + facet_wrap(~aspect) +
  labs(title=TT(paste0("Equidistant group, richness varies: ", paste(EQD, collapse="<"))),
       subtitle="RBH. Recall of truth-enriched terms rises with reference annotation depth",
       x="reference (increasing annotation richness ->)", y="recall (transferred vs truth)",
       fill="evidence set") + theme_bw(base_size=11)
ggsave(file.path(figdir,"fig2_distance_vs_richness.pdf"), p2, width=9, height=4.5)
ggsave(file.path(figdir,"fig2_distance_vs_richness.svg"), p2, width=9, height=4.5)

## ---- Fig 3: IC-bin stratified F1 vs %identity (richness/depth sensitivity) ----
d3 <- m[ic_bin!="all" & list_source=="synthetic" & evset=="all" & method=="rbh",
        .(f1=mean(f1,na.rm=TRUE)), by=.(ref,aspect,ic_bin,median_pident)]
d3[, ic_bin := factor(ic_bin, levels=c("shallow","mid","deep"))]
p3 <- ggplot(d3, aes(median_pident, f1, color=ic_bin, group=ic_bin)) +
  geom_point() + geom_line() + facet_wrap(~aspect) + scale_x_reverse() +
  labs(title=TT("Term specificity (IC) sensitivity (RBH, all-evidence)"),
       subtitle="deep (specific) terms expected to collapse faster than shallow",
       x="median ortholog % identity", y="set-level F1", color="IC bin") + theme_bw(base_size=11)
ggsave(file.path(figdir,"fig3_IC_stratified.pdf"), p3, width=9, height=4)
ggsave(file.path(figdir,"fig3_IC_stratified.svg"), p3, width=9, height=4)

## ---- Semantic similarity (Wang, redundancy-robust) — cleaner main curve ----
s <- fread(file.path(OUT,"metrics/semantic_sim.tsv"))
pid_map <- unique(m[,.(ref,My,median_pident)])
s[, ref := tstrsplit(label,"\\.")[[1]]]
s[, c("method","evset") := tstrsplit(label,"\\.")[2:3]]
s <- merge(s, pid_map, by="ref", all.x=TRUE)
s[, aspect := factor(aspect, levels=c("MF","BP","CC"))]
s[, source := fifelse(grepl("^REAL_", list_id), "real DEG", "synthetic")]
ds <- s[source=="synthetic", .(wang=mean(wang_bma,na.rm=TRUE), se=sd(wang_bma,na.rm=TRUE)/sqrt(.N)),
        by=.(ref,method,evset,aspect,median_pident,My)]
pS <- ggplot(ds[evset=="all"], aes(median_pident, wang, color=aspect, group=aspect)) +
  geom_point(size=2) + geom_line() +
  geom_errorbar(aes(ymin=wang-se, ymax=wang+se), width=1, alpha=.5) +
  facet_wrap(~method, ncol=2) + scale_x_reverse() +
  labs(title=TT("Wang semantic similarity tolerance (all-evidence)"),
       subtitle="redundancy-robust set similarity (trap #5). left=closer. mean +/- SE over synthetic lists",
       x="median ortholog % identity", y="Wang BMA semantic similarity (transferred vs truth set)") +
  theme_bw(base_size=11)
ggsave(file.path(figdir,"fig1_semantic_curve.pdf"), pS, width=9, height=5)
ggsave(file.path(figdir,"fig1_semantic_curve.svg"), pS, width=9, height=5)

## ---- Recall curve (cleaner than F1; richness-sensitive) ----
dr <- m[ic_bin=="all" & list_source=="synthetic" & evset=="all", .(recall=mean(recall,na.rm=TRUE),
        se=sd(recall,na.rm=TRUE)/sqrt(.N)), by=.(method,aspect,median_pident,ref)]
pR <- ggplot(dr, aes(median_pident, recall, color=aspect, group=aspect)) +
  geom_point(size=2)+geom_line()+facet_wrap(~method,ncol=2)+scale_x_reverse()+
  labs(title=TT("Recall vs divergence (all-evidence)"),
       x="median ortholog % identity", y="recall (transferred vs truth)")+theme_bw(base_size=11)
ggsave(file.path(figdir,"fig1c_recall_curve.pdf"), pR, width=9, height=5)
ggsave(file.path(figdir,"fig1c_recall_curve.svg"), pR, width=9, height=5)

## ---- Fig 4: REAL DEG vs SYNTHETIC lists (external validity) — only if real lists exist ----
if (any(s$source=="real DEG")) {
  d4 <- s[method=="rbh" & evset=="all",
          .(wang=mean(wang_bma,na.rm=TRUE), se=sd(wang_bma,na.rm=TRUE)/sqrt(.N), n=.N),
          by=.(source,aspect,median_pident)]
  p4 <- ggplot(d4, aes(median_pident, wang, color=aspect, linetype=source, group=interaction(aspect,source))) +
    geom_point(size=1.8)+geom_line()+scale_x_reverse()+
    labs(title=TT("Real DEGs vs synthetic lists (RBH, all-evidence)"),
         subtitle="real DEGs score low NOT because transfer fails but because EXP-only focal truth is too sparse to validate broad DEG lists (truth sig-set often empty, esp. CC)",
         x="median ortholog % identity", y="Wang semantic similarity", linetype="list source")+
    theme_bw(base_size=11)
  ggsave(file.path(figdir,"fig4_real_vs_synthetic.pdf"), p4, width=9, height=5)
  ggsave(file.path(figdir,"fig4_real_vs_synthetic.svg"), p4, width=9, height=5)
} else {
  cat("[S6] no real DEG lists for this track -> skipping fig4 (synthetic-only)\n")
}

## ---- Guide table: now on Wang semantic similarity (redundancy-robust) ----
gsem <- ds[method=="rbh" & evset=="all", .(ref,My,median_pident,aspect,wang=round(wang,2))]
gf <- m[ic_bin=="all" & list_source=="synthetic" & evset=="all" & method=="rbh",
        .(recall=round(mean(recall,na.rm=TRUE),2)), by=.(ref,aspect)]
g <- merge(gsem, gf, by=c("ref","aspect"))
g[, reliability := fifelse(wang>=0.6,"reliable", fifelse(wang>=0.45,"caution","unreliable"))]
setorder(g, aspect, -median_pident)
fwrite(g, file.path(tabdir,"guide_table.tsv"), sep="\t")

cat("[S6] figures + guide table written\n")
print(dcast(g, ref+My+median_pident~aspect, value.var="wang"))
