# export_tables.R — dump the exact plotted values of every figure as tidy CSVs (for Keynote charts)
suppressMessages({library(dplyr); library(tidyr)})
TAB <- "../data"
OUT <- "../figure_tables"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
rd <- function(f) read.delim(file.path(TAB,f), stringsAsFactors=FALSE, check.names=FALSE)
wr <- function(d,n) write.csv(d, file.path(OUT,paste0(n,".csv")), row.names=FALSE)
ORD <- c("mammal","fish","insect","plant_rice","plant_arabidopsis","fungi")
relab <- c(mammal="Mammal",fish="Fish",insect="Insect",plant_rice="Plant-R",plant_arabidopsis="Plant-A",fungi="Fungi")
lab <- function(x) factor(relab[x], levels=unname(relab[ORD]))   # display label, ordered

sl<-rd("figdata_setlevel.tsv"); ic<-rd("figdata_ic.tsv"); rm0<-rd("figdata_refmeta.tsv")
ci<-rd("crossclade_ci.tsv"); id<-rd("id50.tsv"); cat0<-rd("category.tsv")
wgd<-rd("wgd_fungi.tsv"); strc<-rd("string_conservation.tsv"); exprc<-rd("expr_conservation.tsv")
guide<-rd("guide_table_unified.tsv")
SYN<-sl%>%filter(list_source=="synthetic",evset=="all",size_class=="small")

# Fig1a / S1 / S2 — tolerance curves (Wang + recall, all aspects, with 95% CI)
wr(ci%>%filter(track%in%ORD)%>%transmute(Track=lab(track),aspect,ref,median_pident,
   wang,wang_lo,wang_hi,recall,recall_lo,recall_hi)%>%arrange(aspect,Track,desc(median_pident)),
   "Fig1abc_S1_tolerance_curves")

# Fig1b / S4 — ID50 + floor (+ logistic params)
wr(id%>%filter(track%in%ORD)%>%transmute(Track=lab(track),aspect,ID50=ID50_rep,floor,plateau,
   ID50_logistic,slope_k,fit_R2,identifiable,method_used)%>%arrange(aspect,Track),"Fig1d_S2_ID50_floor")

# Fig2 — distance vs richness (equidistant refs; BP)
perf<-SYN%>%filter(aspect=="BP")%>%group_by(track,ref)%>%summarise(recall=mean(recall,na.rm=TRUE),precision=mean(precision,na.rm=TRUE),.groups="drop")
wr(rm0%>%filter(!is.na(equidistant_group),equidistant_group!="")%>%left_join(perf,by=c("track","ref"))%>%
   transmute(Track=lab(track),ref,equidistant_group,ref_exp_richness,median_pident,recall,precision)%>%arrange(Track,desc(ref_exp_richness)),
   "Fig2_distance_vs_richness")

# Fig3a — conserved vs taxon-specific (BP)
wr(cat0%>%filter(track%in%ORD,aspect=="BP",category%in%c("conserved","specific","other"))%>%
   mutate(wang_bma=as.numeric(wang_bma),median_pident=as.numeric(median_pident))%>%
   group_by(track,category,ref,median_pident)%>%summarise(wang=mean(wang_bma,na.rm=TRUE),.groups="drop")%>%
   transmute(Track=lab(track),category,ref,median_pident,wang)%>%arrange(Track,category,desc(median_pident)),
   "Fig3a_conserved_vs_specific")

# Fig3b / IC — term specificity (BP, RBH, F1+recall)
wr(ic%>%filter(track%in%ORD,method=="rbh",evset=="all",aspect=="BP",ic_bin%in%c("shallow","mid","deep"))%>%
   transmute(Track=lab(track),ic_bin,ref,median_pident=as.numeric(median_pident),f1=as.numeric(f1),recall=as.numeric(recall),n_terms)%>%
   arrange(Track,ic_bin,desc(median_pident)),"Fig3b_IC_stratified")

# Fig4a — best-hit vs RBH (BP)
wr(SYN%>%filter(aspect=="BP",method%in%c("besthit","rbh"))%>%group_by(track,ref,method,median_pident)%>%
   summarise(wang=mean(wang_bma,na.rm=TRUE),.groups="drop")%>%transmute(Track=lab(track),ref,method,median_pident,wang)%>%
   arrange(Track,method,desc(median_pident)),"Fig4a_besthit_vs_rbh")

# Fig4b — RBH penalty (RBH - best-hit) per clade
wr(SYN%>%filter(aspect=="BP",method%in%c("besthit","rbh"))%>%group_by(track,method)%>%
   summarise(across(c(wang_bma,recall,precision),~mean(.,na.rm=TRUE)),.groups="drop")%>%
   pivot_wider(names_from=method,values_from=c(wang_bma,recall,precision))%>%
   transmute(Track=lab(track),Wang_penalty=wang_bma_rbh-wang_bma_besthit,
             precision_penalty=precision_rbh-precision_besthit,recall_penalty=recall_rbh-recall_besthit)%>%arrange(Track),
   "Fig4b_method_penalty")

# Fig4c / S7 — eggNOG: per-ref P-R + mean by method (BP, eval-subset refs)
tri<-SYN%>%group_by(track,ref)%>%filter(n_distinct(method)==3)%>%ungroup()%>%filter(aspect=="BP")
wr(tri%>%group_by(track,ref,method)%>%summarise(recall=mean(recall,na.rm=TRUE),precision=mean(precision,na.rm=TRUE),wang=mean(wang_bma,na.rm=TRUE),.groups="drop")%>%
   transmute(Track=lab(track),ref,method,recall,precision,wang)%>%arrange(Track,ref,method),"Fig4c_S5_eggnog_per_ref")
wr(tri%>%group_by(method)%>%summarise(recall=mean(recall,na.rm=TRUE),precision=mean(precision,na.rm=TRUE),wang=mean(wang_bma,na.rm=TRUE),.groups="drop"),"Fig4c_eggnog_means")

# Fig5a / S5 — orthogonal validation (expression + STRING)
wr(exprc%>%filter(track%in%ORD, ref!="fruitfly")%>%transmute(Track=lab(track),ref,median_pident,expr_conservation,GO_Wang,GO_recall),"Fig5a_orthogonal_expression")  # within-phylum (cross-phylum fruitfly excluded)
wr(strc%>%filter(track%in%ORD)%>%transmute(Track=lab(track),ref,median_pident,string_conservation=suppressWarnings(as.numeric(string_conservation)),GO_Wang,n_genes),"S3_orthogonal_string")

# Fig5b/5c / S6 — robustness
rs<-sl%>%filter(method=="rbh",evset=="all",aspect=="BP")
rs$szbin<-cut(rs$size,breaks=c(0,100,250,500,1000,2500,Inf),labels=c("<100","100-250","250-500","500-1k","1k-2.5k",">2.5k"))
wr(rs%>%filter(!is.na(szbin))%>%group_by(szbin,list_source)%>%summarise(mean_wang=mean(wang_bma,na.rm=TRUE),SE=sd(wang_bma,na.rm=TRUE)/sqrt(sum(!is.na(wang_bma))),n=n(),.groups="drop"),"Fig5b_robustness_sizebin")
tc<-rm0%>%group_by(track)%>%summarise(truth_cov=first(truth_cov),.groups="drop")
wr(rs%>%filter(list_source=="real")%>%group_by(track)%>%summarise(mean_real_wang=mean(wang_bma,na.rm=TRUE),.groups="drop")%>%
   left_join(tc,by="track")%>%transmute(Track=lab(track),truth_cov,mean_real_wang)%>%arrange(truth_cov),"Fig5c_robustness_truthcov")
wr(rs%>%filter(list_source=="synthetic",!is.na(list_ic))%>%transmute(Track=lab(track),size,wang=wang_bma,list_ic),"S4_robustness_size_IC")

# Fig6 / S3 — reliability guide (all aspects, RBH)
wr(guide%>%filter(clade%in%ORD,method=="rbh")%>%transmute(Clade=lab(clade),aspect,distance,wang,n_ref,reliability)%>%arrange(aspect,Clade,distance),"TableS_reliability_full")

# S8 — fungi WGD
wr(wgd,"S6_wgd_fungi")

cat("DONE tables ->",OUT,"\n")
