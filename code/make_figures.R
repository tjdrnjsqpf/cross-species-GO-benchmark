# panels.R — each panel as a SEPARATE file (PDF vector + PNG) for Keynote assembly.
# No patchwork composites. Legends/titles kept minimal; user edits in Keynote.
suppressMessages({library(ggplot2); library(dplyr); library(tidyr)})
TAB <- "../data"
OUT <- "../figures"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
rd <- function(f) read.delim(file.path(TAB,f), stringsAsFactors=FALSE, check.names=FALSE)
savep <- function(p,n,w=5,h=4){ sub<-if(grepl("^Fig",n))"main" else "supp"; d<-file.path(OUT,sub)
  dir.create(d,showWarnings=FALSE,recursive=TRUE)
  p <- p + theme(plot.title=element_blank(), plot.subtitle=element_blank(),
                 axis.title.x=element_text(margin=margin(t=10)),
                 axis.title.y=element_text(margin=margin(r=10)))
  ggsave(file.path(d,paste0(n,".pdf")),p,width=w,height=h) }
ORD <- c("mammal","fish","insect","plant_rice","plant_arabidopsis","fungi")
PAL <- c(mammal="#D95F02",fish="#1B9E77",insect="#B8860B",plant_rice="#E7298A",plant_arabidopsis="#3B6CB7",fungi="#7570B3")
relab <- c(mammal="Mammal",fish="Fish",insect="Insect",plant_rice="Plant-R",plant_arabidopsis="Plant-A",fungi="Fungi")
ORDLAB <- unname(relab[ORD]); lab_track <- as_labeller(relab)
fct <- function(x) factor(x, levels=ORD)
col_t <- scale_color_manual(values=PAL, breaks=ORD, labels=ORDLAB)
fil_t <- scale_fill_manual(values=PAL, breaks=ORD, labels=ORDLAB, guide="none")

sl  <- rd("figdata_setlevel.tsv")
ic  <- rd("figdata_ic.tsv")
rm0 <- rd("figdata_refmeta.tsv")
ci  <- rd("crossclade_ci.tsv") %>% filter(track%in%ORD) %>% mutate(track=fct(track))
id  <- rd("id50.tsv") ; id$identifiable <- id$identifiable %in% c("True","TRUE",TRUE); id <- id%>%filter(track%in%ORD)%>%mutate(track=fct(track))
cat0<- rd("category.tsv")
wgd <- rd("wgd_fungi.tsv")
strc<- rd("string_conservation.tsv"); strc$string_conservation<-suppressWarnings(as.numeric(strc$string_conservation))
exprc<-rd("expr_conservation.tsv")
SYN <- sl %>% filter(list_source=="synthetic", evset=="all", size_class=="small")

# ---------- tolerance curves: one panel per aspect (Wang + recall) ----------
tol <- function(asp, ycol, ylo, yhi, ylab, tag){
  d <- ci %>% filter(aspect==asp)
  ggplot(d, aes(median_pident, .data[[ycol]], color=track, fill=track)) +
    geom_ribbon(aes(ymin=.data[[ylo]], ymax=.data[[yhi]]), alpha=.12, color=NA) +
    geom_line(alpha=.85)+geom_point(size=1.2,alpha=.9)+col_t+fil_t+scale_x_reverse()+coord_cartesian(ylim=c(0,1))+
    labs(title=paste0(tag," (",asp,")"), x="median ortholog % identity", y=ylab, color="track")+
    theme_bw(11)+theme(legend.position="right")
}
wnm <- c(BP="Fig1a_tolerance_wang_BP", MF="Fig1b_tolerance_wang_MF", CC="Fig1c_tolerance_wang_CC")
rnm <- c(BP="S1a_tolerance_recall_BP", MF="S1b_tolerance_recall_MF", CC="S1c_tolerance_recall_CC")
for(a in c("BP","MF","CC")){
  savep(tol(a,"wang","wang_lo","wang_hi","Wang semantic similarity","Tolerance"), wnm[[a]],6,4)
  savep(tol(a,"recall","recall_lo","recall_hi","set-level recall","Tolerance recall"), rnm[[a]],6,4)
}

# ---------- ID50 + floor bars ----------
idb <- id%>%filter(aspect=="BP")%>%transmute(track=fct(track), `ID50 (% identity)`=ID50_rep, `floor (Wang)`=floor)%>%
  pivot_longer(-track,names_to="metric",values_to="v")%>%mutate(metric=factor(metric,levels=c("ID50 (% identity)","floor (Wang)")))
savep(ggplot(idb,aes(track,v,fill=track))+geom_col()+facet_wrap(~metric,scales="free_y")+fil_t+
  scale_x_discrete(labels=relab)+labs(x=NULL,y=NULL)+
  theme_bw(11)+theme(axis.text.x=element_text(angle=40,hjust=1),legend.position="none"), "Fig1d_ID50_floor_BP",7,3.8)

# ---------- ID50 logistic fits (BP) ----------
cb<-ci%>%filter(aspect=="BP"); idbb<-id%>%filter(aspect=="BP"); xr<-range(cb$median_pident,na.rm=TRUE)
crv<-idbb%>%rowwise()%>%do({r<-.;x<-seq(xr[1],xr[2],length.out=200)
  data.frame(track=r$track,x=x,y=r$floor+(r$plateau-r$floor)/(1+exp(-r$slope_k*(x-r$ID50_logistic))),identifiable=r$identifiable)})%>%ungroup()%>%mutate(track=fct(track))
savep(ggplot()+geom_point(data=cb,aes(median_pident,wang,color=track),size=1.5,alpha=.5)+
  geom_line(data=crv,aes(x,y,color=track,linetype=identifiable),linewidth=1)+
  geom_vline(data=idbb%>%filter(identifiable),aes(xintercept=ID50_rep,color=track),linetype=3,linewidth=.5,show.legend=FALSE)+
  col_t+scale_linetype_manual(values=c(`TRUE`="solid",`FALSE`="22"),labels=c(`TRUE`="identifiable",`FALSE`="NE->interp"),name="fit")+
  scale_x_reverse()+labs(title="ID50 logistic fits (BP)",x="median ortholog % identity",y="Wang (BP)",color="track")+theme_bw(11),
  "S2_ID50_logistic_BP",7,4.5)

# ---------- reliability guide -> now Table 1 (no heatmap figure); see export_tables.R ----------

# ---------- distance vs richness: recall & precision (separate) ----------
eqd<-rm0%>%filter(!is.na(equidistant_group),equidistant_group!="")%>%mutate(track=fct(track))
perf<-SYN%>%filter(aspect=="BP")%>%group_by(track,ref)%>%summarise(recall=mean(recall,na.rm=TRUE),precision=mean(precision,na.rm=TRUE),.groups="drop")
dr<-eqd%>%left_join(perf,by=c("track","ref"))
drp<-function(ycol,ylab,xx,hj){
  rr<-suppressWarnings(cor.test(dr$ref_exp_richness, dr[[ycol]], method="spearman"))
  ann<-sprintf("Spearman rho = %+.2f, p = %.3f", as.numeric(rr$estimate), rr$p.value)
  ggplot(dr,aes(ref_exp_richness,.data[[ycol]],color=track))+geom_point(size=2.3,alpha=.85)+
  geom_smooth(method="lm",se=FALSE,linewidth=.5,aes(group=1),color="grey40",linetype=2)+col_t+
  annotate("text",x=xx,y=Inf,label=ann,hjust=hj,vjust=1.6,size=3.3)+
  labs(x="reference richness (EXP GO terms/gene)",y=paste0(ylab," (BP)"),color="track")+theme_bw(11) }
savep(drp("recall","recall",-Inf,-0.05),"Fig2a_distance_richness_recall",5.5,4)      # label top-left (empty)
savep(drp("precision","precision",Inf,1.05),"Fig2b_distance_richness_precision",5.5,4) # label top-right (empty)

# ---------- conserved vs specific (track-faceted) ----------
cc<-cat0%>%filter(track%in%ORD,aspect=="BP",category%in%c("conserved","specific","other"))%>%
  mutate(wang_bma=as.numeric(wang_bma),median_pident=as.numeric(median_pident))%>%
  group_by(track,category,ref,median_pident)%>%summarise(wang=mean(wang_bma,na.rm=TRUE),.groups="drop")%>%
  mutate(track=fct(track),category=factor(category,levels=c("conserved","specific","other")))
savep(ggplot(cc,aes(median_pident,wang,color=category))+geom_point(size=1.1,alpha=.6)+
  geom_smooth(method="loess",se=FALSE,span=1,linewidth=.8)+facet_wrap(~track,nrow=2,labeller=lab_track)+
  scale_color_manual(values=c(conserved="#2166AC",specific="#B2182B",other="grey60"))+scale_x_reverse()+coord_cartesian(ylim=c(0,.95))+
  labs(title="Conserved vs taxon-specific (BP, RBH)",x="median ortholog % identity",y="Wang (BP)",color="GO category")+
  theme_bw(11)+theme(legend.position="bottom"),"Fig3a_conserved_vs_specific",9,5.2)

# ---------- IC stratified (track-faceted) ----------
icd<-ic%>%filter(track%in%ORD,method=="rbh",evset=="all",aspect=="BP",ic_bin%in%c("shallow","mid","deep"))%>%
  mutate(track=fct(track),ic_bin=factor(ic_bin,levels=c("shallow","mid","deep")),f1=as.numeric(f1),median_pident=as.numeric(median_pident))
savep(ggplot(icd,aes(median_pident,f1,color=ic_bin))+geom_point(size=.9,alpha=.45)+
  geom_smooth(method="loess",se=FALSE,span=1,linewidth=.8)+facet_wrap(~track,nrow=2,labeller=lab_track)+scale_x_reverse()+coord_cartesian(ylim=c(0,NA))+
  scale_color_manual(values=c(shallow="#4575B4",mid="#F4A100",deep="#D73027"))+
  labs(title="Term specificity (IC) sensitivity (BP, RBH, F1)",x="median ortholog % identity",y="F1 (BP)",color="IC bin")+
  theme_bw(11)+theme(legend.position="bottom"),"Fig3b_IC_stratified",9,5.2)

# ---------- besthit vs rbh (method-faceted) ----------
bv<-SYN%>%filter(aspect=="BP",method%in%c("besthit","rbh"))%>%group_by(track,ref,method,median_pident)%>%
  summarise(wang=mean(wang_bma,na.rm=TRUE),.groups="drop")%>%mutate(track=fct(track))
savep(ggplot(bv,aes(median_pident,wang,color=track))+geom_line(alpha=.7)+geom_point(size=1,alpha=.85)+
  facet_wrap(~method,nrow=1)+col_t+scale_x_reverse()+labs(title="best-hit vs RBH (BP)",x="median ortholog % identity",y="Wang (BP)",color="track")+
  theme_bw(11),"Fig4a_besthit_vs_rbh",9,4)

# ---------- method penalty (metric-faceted) ----------
mp<-SYN%>%filter(aspect=="BP",method%in%c("besthit","rbh"))%>%group_by(track,method)%>%
  summarise(across(c(wang_bma,recall,precision),~mean(.,na.rm=TRUE)),.groups="drop")%>%
  pivot_wider(names_from=method,values_from=c(wang_bma,recall,precision))%>%
  transmute(track=fct(track),Wang=wang_bma_rbh-wang_bma_besthit,precision=precision_rbh-precision_besthit,recall=recall_rbh-recall_besthit)%>%
  pivot_longer(-track,names_to="metric",values_to="penalty")%>%mutate(metric=factor(metric,levels=c("Wang","precision","recall")))
savep(ggplot(mp,aes(track,penalty,fill=track))+geom_col()+facet_wrap(~metric,nrow=1,scales="free_y")+fil_t+scale_x_discrete(labels=relab)+geom_hline(yintercept=0,linewidth=.3)+
  labs(title="RBH penalty (RBH - best-hit; negative = RBH worse)",x=NULL,y="RBH - best-hit")+
  theme_bw(11)+theme(axis.text.x=element_text(angle=40,hjust=1)),"Fig4b_method_penalty",10,4)

# ---------- eggnog: PR scatter + bars (separate) ----------
tri<-SYN%>%group_by(track,ref)%>%filter(n_distinct(method)==3)%>%ungroup()%>%filter(aspect=="BP")
trr<-tri%>%group_by(track,ref,method)%>%summarise(recall=mean(recall,na.rm=TRUE),precision=mean(precision,na.rm=TRUE),.groups="drop")
savep(ggplot(trr,aes(recall,precision,color=method))+geom_line(aes(group=interaction(track,ref)),color="grey80",linewidth=.3)+
  geom_point(size=2,alpha=.85)+theme_bw(11)+labs(title="P-R per reference (BP)",x="recall",y="precision"),"S5_eggnog_PR",5.5,4)
trw<-trr%>%left_join(tri%>%group_by(track,ref,method)%>%summarise(wang=mean(wang_bma,na.rm=TRUE),.groups="drop"),by=c("track","ref","method"))%>%
  group_by(method)%>%summarise(across(c(recall,precision,wang),list(m=~mean(.,na.rm=TRUE),se=~sd(.,na.rm=TRUE)/sqrt(sum(!is.na(.)))),.names="{.col}_{.fn}"),.groups="drop")%>%
  pivot_longer(-method,names_to=c("metric",".value"),names_sep="_")
savep(ggplot(trw,aes(metric,m,fill=method))+geom_col(position=position_dodge(.8),width=.7)+
  geom_errorbar(aes(ymin=m-se,ymax=m+se),position=position_dodge(.8),width=.2)+
  theme_bw(11)+labs(x=NULL,y="mean (BP)  +/- SE across refs"),"Fig4c_eggnog_bars",5.5,4)

# ---------- wgd fungi (metric-faceted) ----------
wl<-wgd%>%group_by(wgd)%>%summarise(frac_nonrecip=mean(frac_nonrecip,na.rm=TRUE),many_to_one=mean(many_to_one,na.rm=TRUE),.groups="drop")%>%
  filter(wgd%in%c("post","pre"))%>%pivot_longer(c(frac_nonrecip,many_to_one),names_to="metric",values_to="v")
savep(ggplot(wl,aes(wgd,v,fill=wgd))+geom_col()+facet_wrap(~metric,nrow=1)+scale_fill_manual(values=c(post="#D95F02",pre="#1B9E77"),guide="none")+
  labs(title="Fungi post- vs pre-WGD paralog ambiguity (raw; distance NOT controlled)",x=NULL,y="fraction")+theme_bw(11),"S6_wgd_fungi",6,4)

# ---------- S17 expression (track-faceted) + string ----------
ex2<-exprc%>%filter(!is.na(expr_conservation),!is.na(GO_Wang),ref!="fruitfly")%>%group_by(track)%>%filter(n()>=3)%>%ungroup()%>%mutate(track=fct(track))  # within-phylum only: drop cross-phylum fruitfly ref + insect(n=1,NE)
st2<-strc%>%filter(!is.na(string_conservation),!is.na(GO_Wang),n_genes>=50)%>%mutate(track=fct(track))
rho<-function(d,xc) d%>%group_by(track)%>%summarise(txt=sprintf("rho=%+.2f (n=%d)",suppressWarnings(cor(.data[[xc]],GO_Wang,method="spearman")),n()),.groups="drop")
savep(ggplot(ex2,aes(expr_conservation,GO_Wang,color=track))+geom_point(size=2.3,alpha=.85,show.legend=FALSE)+
  geom_smooth(method="lm",se=FALSE,color="grey40",linetype=2,linewidth=.5)+col_t+
  geom_text(data=rho(ex2,"expr_conservation"),aes(x=Inf,y=Inf,label=txt),hjust=1.05,vjust=1.4,size=3,color="black")+
  facet_wrap(~track,scales="free",nrow=1,labeller=lab_track)+labs(title="Bgee expression conservation vs GO-transfer (animals)",x="expression conservation",y="GO-transfer Wang (BP)")+theme_bw(11),"Fig5a_S17_expression",10,3.4)
savep(ggplot(st2,aes(string_conservation,GO_Wang,color=track))+geom_point(size=1.8,alpha=.8,show.legend=FALSE)+
  geom_smooth(method="lm",se=FALSE,color="grey40",linetype=2,linewidth=.5)+col_t+
  geom_text(data=rho(st2,"string_conservation"),aes(x=Inf,y=Inf,label=txt),hjust=1.05,vjust=1.4,size=3,color="black")+
  facet_wrap(~track,scales="free",nrow=2,labeller=lab_track)+labs(title="STRING neighbour Jaccard vs GO-transfer (inconclusive)",x="STRING neighbour Jaccard",y="GO-transfer Wang (BP)")+theme_bw(11),"S3_S17_string",9,5.5)

# ---------- robustness D1/D2/D3 (separate) ----------
rs<-sl%>%filter(method=="rbh",evset=="all",aspect=="BP")%>%mutate(list_source=factor(list_source))
rs$szbin<-cut(rs$size,breaks=c(0,100,250,500,1000,2500,Inf),labels=c("<100","100-250","250-500","500-1k","1k-2.5k",">2.5k"))
binm<-rs%>%filter(!is.na(szbin))%>%group_by(szbin,list_source)%>%summarise(wang=mean(wang_bma,na.rm=TRUE),se=sd(wang_bma,na.rm=TRUE)/sqrt(sum(!is.na(wang_bma))),n=n(),.groups="drop")
savep(ggplot(binm,aes(szbin,wang,fill=list_source))+geom_col(position=position_dodge(.8),width=.7)+
  geom_errorbar(aes(ymin=wang-se,ymax=wang+se),position=position_dodge(.8),width=.2)+
  labs(title="Wang by size bin (real < synthetic at matched size)",x="list size bin (#genes)",y="mean Wang (BP) +/- SE",fill="list source")+
  theme_bw(11)+theme(axis.text.x=element_text(angle=25,hjust=1)),"Fig5b_robustness_D1_sizebin",6,4)
tc<-rm0%>%group_by(track)%>%summarise(truth_cov=first(truth_cov),.groups="drop")
rw<-rs%>%filter(list_source=="real")%>%group_by(track)%>%summarise(wang_real=mean(wang_bma,na.rm=TRUE),.groups="drop")%>%left_join(tc,by="track")%>%mutate(track=fct(track))
rwr<-suppressWarnings(cor.test(rw$truth_cov,rw$wang_real,method="spearman"))
savep(ggplot(rw,aes(truth_cov,wang_real,color=track))+geom_point(size=3)+col_t+
  ggrepel::geom_text_repel(aes(label=relab[as.character(track)]),size=3,show.legend=FALSE)+
  annotate("text",x=-Inf,y=Inf,hjust=-0.08,vjust=1.6,size=3.4,
           label=sprintf("Spearman rho = %+.2f, p = %.3f", as.numeric(rwr$estimate), rwr$p.value))+
  labs(title="real-DEG Wang vs truth coverage",x="focal truth coverage",y="mean real-DEG Wang (BP)")+theme_bw(11)+theme(legend.position="none"),"Fig5c_robustness_D2_truthcov",5.5,4)
synic<-rs%>%filter(list_source=="synthetic",!is.na(list_ic))
savep(ggplot(synic,aes(size,wang_bma,color=list_ic))+geom_point(size=.8,alpha=.5)+scale_x_log10()+scale_color_gradient(low="#FDE725",high="#440154")+
  labs(title="synthetic Wang vs size, by IC",x="list size (#genes, log)",y="Wang (BP)",color="list IC")+theme_bw(11),"S4_robustness_D3_sizeIC",5.5,4)

cat("DONE panels ->", OUT, "\n")
