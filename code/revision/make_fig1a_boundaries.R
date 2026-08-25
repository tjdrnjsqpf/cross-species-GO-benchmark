# Revision R1-9 — Fig 1A tolerance panels with clade-boundary markers.
# Reproduces make_figures.R Fig1a/b/c aesthetics exactly, adding one dashed
# vertical line per focal species at the identity of its phylogenetically
# nearest OUT-of-clade reference ("first outgroup" = clade boundary).
# NOTE (legend text): for fish and insect one within-clade reference (grass
# carp; pea aphid) lies beyond the first outgroup on the identity axis —
# taxonomy and sequence identity are not perfectly collinear, which is why
# the divergence axis uses identity.
suppressMessages({library(ggplot2); library(dplyr)})
# resolve repo root from this script location (code/revision/)
args <- commandArgs(trailingOnly=FALSE)
sp <- sub("^--file=", "", args[grep("^--file=", args)])
BASE <- normalizePath(file.path(dirname(sp), "../.."))
TAB  <- file.path(BASE, "data")
OUT  <- file.path(BASE, "data/revision")
TABD <- file.path(BASE, "data/revision")
rd <- function(f) read.delim(file.path(TAB,f), stringsAsFactors=FALSE, check.names=FALSE)

ORD <- c("mammal","fish","insect","plant_rice","plant_arabidopsis","fungi")
PAL <- c(mammal="#D95F02",fish="#1B9E77",insect="#B8860B",
         plant_rice="#E7298A",plant_arabidopsis="#3B6CB7",fungi="#7570B3")
relab <- c(mammal="Mammal",fish="Fish",insect="Insect",
           plant_rice="Plant-R",plant_arabidopsis="Plant-A",fungi="Fungi")
fct <- function(x) factor(x, levels=ORD)
col_t <- scale_color_manual(values=PAL, breaks=ORD, labels=unname(relab[ORD]))
fil_t <- scale_fill_manual(values=PAL, breaks=ORD, labels=unname(relab[ORD]), guide="none")

ci <- rd("crossclade_ci.tsv") %>% filter(track %in% ORD) %>% mutate(track=fct(track))

# clade boundary = identity of the nearest out-of-clade reference (fig1a_clade_boundaries.tsv)
bnd <- read.delim(file.path(TABD, "fig1a_clade_boundaries.tsv"), stringsAsFactors=FALSE) %>%
  transmute(track=fct(track), xb=pid_out)

tol <- function(asp){
  d <- ci %>% filter(aspect==asp)
  ggplot(d, aes(median_pident, wang, color=track, fill=track)) +
    geom_vline(data=bnd, aes(xintercept=xb, color=track),
               linetype="22", linewidth=.45, alpha=.75, show.legend=FALSE) +
    geom_ribbon(aes(ymin=wang_lo, ymax=wang_hi), alpha=.12, color=NA) +
    geom_line(alpha=.85) + geom_point(size=1.2, alpha=.9) +
    col_t + fil_t + scale_x_reverse() + coord_cartesian(ylim=c(0,1)) +
    labs(x="median ortholog % identity", y="Wang semantic similarity", color="track") +
    theme_bw(11) + theme(legend.position="right",
                         axis.title.x=element_text(margin=margin(t=10)),
                         axis.title.y=element_text(margin=margin(r=10)))
}
nm <- c(BP="Fig1a_tolerance_wang_BP_boundaries",
        MF="Fig1b_tolerance_wang_MF_boundaries",
        CC="Fig1c_tolerance_wang_CC_boundaries")
for(a in names(nm)){
  p <- tol(a)
  ggsave(file.path(OUT, paste0(nm[[a]], ".pdf")), p, width=6, height=4)
  ggsave(file.path(OUT, paste0(nm[[a]], ".png")), p, width=6, height=4, dpi=300)
}
cat("wrote", paste(nm, collapse=", "), "\n")
