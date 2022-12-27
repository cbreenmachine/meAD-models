
suppressPackageStartupMessages({
    library(viridis)
    library(tidyverse)
    library(data.table)
    library(argparse)
    library(cowplot)
    library(ggsci)
    library(ggbio)
    library(EnsDb.Hsapiens.v86)
    library(clusterProfiler)
}) 

# Gets us plot_GO()
source("../functions.R")

#--> Ensembl database and filtered for coding genes only
ensdb <- EnsDb.Hsapiens.v86
seqlevelsStyle(ensdb) <- "UCSC"

parser <- ArgumentParser(description='')
parser$add_argument('--ifile', default = "../../dataDerived/analysis-controlLOAD/test-diagnostic-group-coded/experimentSummary/pvals.bed", help = "input BED file")
parser$add_argument('--odir', default = "../../figs/2022-paper/dmGenes/", help = "input BED file")
parser$add_argument('--upstream', default = 5000, help = "number of bases upstream from TSS")
parser$add_argument('--downstream', default = 300, help = "number of bases downstream from transcription stop")
parser$add_argument('--effect_cut', default = 0.5, help = "number of bases downstream from transcription stop")
parser$add_argument('--lfdr_cut', default = 0.05, help = "number of bases downstream from transcription stop")
parser$add_argument('--gene_p_cut', default = 0.01, help = "number of bases downstream from transcription stop")
args <- parser$parse_args()

# Derived constants
UP <- as.numeric(args$upstream)
DOWN <- as.numeric(args$downstream)
EFFECT.CUT <- as.numeric(args$effect_cut)
LFDR.CUT <- as.numeric(args$lfdr_cut)
GENE.P.CUT <- as.numeric(args$gene_p_cut)

odir <- file.path(args$odir, paste0("lfdr", LFDR.CUT, "-effect", EFFECT.CUT, "-upstream", UP))
dir.create(odir, showWarn=F, recurs=T)

# Something like dmps.all
# prefix <- str_remove(basename(args$ifile), ".bed")

# Read data and convert to GRanges
data <- fread(args$ifile)
dmps <- dplyr::filter(data, lfdr < LFDR.CUT, abs(pi.diff) > EFFECT.CUT)

# Convert to GRanges for overlap
data.gr <- makeGRangesFromDataFrame(data, keep.extra.columns=T, starts.in.df.are.0based = T)
dmps.gr <- makeGRangesFromDataFrame(dmps, keep.extra.columns=T, starts.in.df.are.0based = T)

seqlevelsStyle(data.gr) <- "UCSC"
seqlevelsStyle(dmps.gr) <- "UCSC"


######################################
######### ANNOTATE TO GENES ##########
######################################

#TODO: fix out of range warning
expand_genes <- function(gr, upstream, downstream) {
# https://bioinformatics.stackexchange.com/questions/4390/expand-granges-object-different-amounts-upstream-vs-downstream
  strand_is_minus = strand(gr) == "-"
  on_plus = which(!strand_is_minus)
  on_minus = which(strand_is_minus)
  start(gr)[on_plus] = start(gr)[on_plus] - upstream
  start(gr)[on_minus] = start(gr)[on_minus] - downstream
  end(gr)[on_plus] = end(gr)[on_plus] + downstream
  end(gr)[on_minus] = end(gr)[on_minus] + upstream
  return(gr)
}

# Filter and expand ensdb
ensdb.subset <- genes(ensdb, filter = GeneBiotypeFilter('protein_coding'))
ensdb.expanded <- expand_genes(ensdb.subset, upstream=UP, downstream=DOWN)

# Overlap our loci and ENSDB
all.overlaps <- findOverlaps(data.gr, ensdb.expanded)
dmp.overlaps <- findOverlaps(dmps.gr, ensdb.expanded)

overlaps_to_counts <- function(overlaps){
  # Subset ensdb for genes with sig CpGs (DMPs)
  # There will be multiple duplicates because 
  # e.g. 10 DMPs will all overlap the same gene
  tmp <- ensdb.expanded[subjectHits(overlaps)]

  # How to handle One gene_name mapping to multiple ENSEMBL IDs?
  # Tally up the number of CpGs (N or X) overlapping gene
  overlaps.df <- data.frame(mcols(tmp)) %>% 
    group_by(gene_name) %>%
    summarize(count = n())
  overlaps.df
}

# Each gene represents a binomial experiment
# with N = (# CpGs in gene) trials and
# X = (# DMPs) successes. This is realized, so we use 'x'
N.df <- overlaps_to_counts(all.overlaps) %>% rename("N" = "count")
x.df <- overlaps_to_counts(dmp.overlaps) %>% rename("x" = "count")
n.genes.tested <- nrow(N.df)

# By dropping NAs, we're implicitly getting rid of genes
# where there ARE CpGs, but zero DMPs. We won't 
# calculate p-values here for speed, but will correct for them in 
# Bonferroni (i.e. divide alpha by ~20k and not just 500)
test.df <- full_join(x.df, N.df, by="gene_name") %>% drop_na()

# Empirical Bayes style 
global.p <- nrow(dmps) / nrow(data)

# Test one-sided for genes who have more significant DMPs than expected
test.df$p.binom <- NA
for (i in 1:nrow(test.df)){
  # This can be vectorizes, but out of 
  # abundance of caution...
  N <- test.df$N[i]
  x <- test.df$x[i]

  # Run test and add 
  test.out <- binom.test(x=x, n=N, p=global.p, alternative = "greater")
  test.df$p.binom[i] <- test.out[[3]]
}

# Which genes are significant
sig.ix <- (test.df$p.binom < (GENE.P.CUT / n.genes.tested))
test.df[sig.ix, ] 

name.ix <- (ensdb.expanded$gene_name %in% test.df$gene_name)
names.df <- data.frame(mcols(ensdb.subset[name.ix, c("gene_name", "gene_id")]))
test.id.df <- left_join(test.df, names.df, by = "gene_name") %>% 
  arrange(gene_name) %>% 
  dplyr::filter(!duplicated(gene_name))

ofile <- file.path(odir, "genesTested.csv")
write_csv(test.id.df, ofile)


# Just significant genes
out.df <- test.id.df  %>%
  dplyr::filter(p.binom < (GENE.P.CUT / n.genes.tested))

ofile <- file.path(odir, "genesTestedSig.txt")
write.table(sort(out.df$gene_name), ofile,row=F, quote=F, col.names=F)


########################################
########## GENE ONTOLOGY ###############
########################################

go.out <- clusterProfiler::enrichGO(out.df$gene_id, OrgDb = 'org.Hs.eg.db', 
    keyType = "ENSEMBL", ont="ALL", pAdjustMethod="fdr")

p <- plot_GO(go.out)

ofile <- file.path(odir,  "DMGenesGO.png")
cowplot::save_plot(ofile, p, base_width=15.5, base_height=9)
