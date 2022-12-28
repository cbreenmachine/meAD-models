# 3-exportDMPsAsBED.R
suppressPackageStartupMessages({
    library(data.table)
    library(argparse)
    library(DSS)
    library(parallel)
    library(GenomicRanges)
    library(EnsDb.Hsapiens.v86)
})

source("ensembl-genes.functions.R")

ensdb.expanded <- get_expanded_ensdb(5000, 300)

mean.effect.thresh <- 0.1
max.effect.thresh <- 0.5
lfdr.cut <- 0.05


#TODO: 
#Additional filter--at least one CpG in upper 75% of effect size
#Or average
parser <- ArgumentParser()
parser$add_argument("--ifile", default= "../DataDerived/ControlLOAD/test-diagnostic-group-coded/ExperimentSummary/pvals.bed", help='Where are the models stored')
parser$add_argument("--odir", default= "../DataDerived/ControlLOAD/test-diagnostic-group-coded/ExperimentSummary/", help='Where are the models stored')
args <- parser$parse_args()

dir.create(args$odir, recursive=T)

ofile.1 <- file.path(args$odir, "DMRegions.bed")
ofile.2 <- file.path(args$odir, "DMRegions-filtered.bed")

# Read data and get into DSS-friendly form
df <- fread(args$ifile) 
df$pos <- df$start
df$pvals <- df$p.corrected
df$fdrs <- df$lfdr

# Call DMRs
p.thresh <- 0.01 / nrow(df) # default, could make this smaller
pct.sig <- 0.99 # default is 0.5, so this is conservative

# May not anchor ro summary
global.effects <- abs(df$pi.diff)
summary(global.effects)

# Get the columns in the correct (DSS) order
input.df <- dplyr::select(df, c("chr", "pos", "stat", "pvals", "lfdr"))
class(input.df)[2] <- "DMLtest.multiFactor"

# Use DSS's DMR findeer
dmrs <- DSS::callDMR(input.df, pct.sig = pct.sig)

add_more_info_to_dmrs <- function(i){
    # Pull out the needed DMR ifnormation
    dmr.chr <- dmrs$chr[i]
    dmr.start <- dmrs$start[i]
    dmr.end <- dmrs$end[i]
    
    cpgs <- df %>%
        dplyr::filter(chr == dmr.chr, start >= dmr.start, start <= dmr.end)

    is.persistent <- (median(abs(cpgs$pi.diff)) > mean.effect.thresh)
    has.one.large <- (max(abs(cpgs$pi.diff)) > max.effect.thresh)
    n.sig <- sum(cpgs$lfdr < lfdr.cut)
    n.sig.hyper <- sum(cpgs$lfdr < lfdr.cut & cpgs$pi.diff > 0)

    return(c(is.persistent, has.one.large, n.sig, n.sig.hyper))
}

# Takes a minute to compute this
out <- mclapply(FUN=add_more_info_to_dmrs, 1:nrow(dmrs), mc.cores=12)
filter.flags <- data.frame(do.call(rbind, out)) %>% 
    dplyr::rename("persistent.effect" = "X1", "one.large.effect" = "X2", "n.sig" = "X3", "n.sig.hyper" = "X4")

# All DMRs get written
dmrs.all <- cbind(dmrs, filter.flags) 

dmrs.gr <- makeGRangesFromDataFrame(dmrs.all, keep.extra.columns=T)
genes.ix <- nearest(dmrs.gr, ensdb.expanded, ignore.strand=T)
dd <- distanceToNearest(dmrs.gr, ensdb.expanded, ignore.strand = T)

# Add nearest gene_name and distane to dmrs dataframe
dmrs.all$nearest_gene_name <- ensdb.expanded[genes.ix]$gene_name
dmrs.all$nearest_gene_id <- ensdb.expanded[genes.ix]$gene_id
dmrs.all$dist_to_nearest_gene <- data.frame(dd)$distance

# Nearest gene and distance TO/FROM nearest gene
fwrite(dmrs.all, ofile.1, sep="\t")

# And the ones that have strong effect size
dmrs.filt <- dmrs.all %>%
    dplyr::filter(persistent.effect | one.large.effect)

fwrite(dmrs.filt, ofile.2, sep= "\t")
#END