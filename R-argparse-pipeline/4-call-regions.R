# 3-exportDMPsAsBED.R
suppressPackageStartupMessages({
    library(data.table)
    library(argparse)
    library(magrittr)
    library(DSS)
    library(parallel)
    library(GenomicRanges)
    library(EnsDb.Hsapiens.v86)
})

source("ensembl-genes.functions.R")

parser <- ArgumentParser()
parser$add_argument("--ifile", default= "../DataDerived/ControlLOAD/test-diagnostic-group-coded/Summaries/pvals.bed", help='Where are the models stored')
parser$add_argument("--odir", default= "../DataDerived/ControlLOAD/test-diagnostic-group-coded/Summaries/")
parser$add_argument("--upstream", default = 5000)
parser$add_argument("--downstream", default = 300)

parser$add_argument("--lfdr", default = 0.05)
# parser$add_argument("--percent_significant", default = 0.5, help = "Passed to DSS callDMR() function")
# parser$add_argument("--min_cpgs", default = 3, help = "Passed to DSS callDMR() function")
# Migrate to running in R, so that we are not loading data over and over again.

args <- parser$parse_args()

PERCENT_SIG_RANGE <- c(0.5, 0.6, 0.7, 0.8, 0.9)
MIN_CG_RANGE <- c(3, 5, 7, 9, 15)


######################################
###### Command line arguments ########
######################################

# Set the upstream and downstream boundaries based on 
# command line arguments
ensdb.expanded <- suppressWarnings(get_expanded_ensdb(args$upstream, args$downstream))
lfdr.cut <- args$lfdr

# ODIR is where DMRegions.bed will be written
# all the deviations from the default parameters will
# be written in the SWEEP.ODIR
ODIR <- args$odir
SWEEP.ODIR <- file.path(ODIR, "DMRsParamSweep")
dir.create(SWEEP.ODIR, recursive=T, showWarn=F)


# Read data and get into DSS-friendly form
df <- fread(args$ifile) 
df$pos <- df$start
df$pvals <- df$p.from.DSS
df$lfdr <- df$lfdr.from.ss # shouldnt' effect the outcome

# May not anchor ro summary
# global.effects <- abs(df$pi.diff)

# Get the columns in the correct (DSS) order
input.df <- dplyr::select(df, c("chr", "pos", "stat", "pvals", "lfdr"))
class(input.df)[2] <- "DMLtest.multiFactor"


# if (is.null(dmrs)){
#     file.create(ofile)
#     stop("No DMRs")
# }

summarize_one_dmr <- function(i){
    # Pull out the needed DMR ifnormation
    dmr.chr <- dmrs$chr[i]
    dmr.start <- dmrs$start[i]
    dmr.end <- dmrs$end[i]
    
    cpgs <- df %>%
        dplyr::filter(chr == dmr.chr, start >= dmr.start, start <= dmr.end)

    median.effect <- median(abs(cpgs$pi.diff))
    max.effect <- max(abs(cpgs$pi.diff))
    n.sig <- sum(cpgs$lfdr < lfdr.cut)
    n.sig.hyper <- sum(cpgs$lfdr < lfdr.cut & cpgs$pi.diff > 0)

    return(c(median.effect, max.effect, n.sig, n.sig.hyper))
}

# Takes a minute to compute this
summarize_all_dmrs <- function(dmrs){
    # Call the 
    info <- mclapply(FUN=summarize_one_dmr, 1:nrow(dmrs), mc.cores=12)
    flags <- data.frame(do.call(rbind, info)) %>% 
        dplyr::rename("persistent.effect" = "X1", 
                    "one.large.effect" = "X2", 
                    "n.sig" = "X3", 
                    "n.sig.hyper" = "X4")

        
    # All DMRs get written
    cbind(dmrs, flags) 
}



add_gene_information <- function(dmrs){
    # Cast to GRanges for ease of use; get closest gene
    dmrs.gr <- makeGRangesFromDataFrame(dmrs, keep.extra.columns=T)
    genes.ix <- nearest(dmrs.gr, ensdb.expanded, ignore.strand=T)
    dd <- distanceToNearest(dmrs.gr, ensdb.expanded, ignore.strand = T)

    # Add nearest gene_name and distane to dmrs dataframe
    dmrs$nearest_gene_name <- ensdb.expanded[genes.ix]$gene_name
    dmrs$nearest_gene_id <- ensdb.expanded[genes.ix]$gene_id
    dmrs$dist_to_nearest_gene <- data.frame(dd)$distance
    dmrs
}

make_oname <- function(pct.sig, min.cg){
    pct.sig.string <- pct.sig * 100
    zz <- paste0("DMRegions-minCG-", min.cg, "-pctSig-", pct.sig.string, ".bed")
    file.path(SWEEP.ODIR, zz)
}

# Store the numnber of DMRs for each parameter configuration
N.configs <- length(PERCENT_SIG_RANGE) * length(MIN_CG_RANGE)

summary.df <- data.frame(
    PercentSig = rep(NA, N.configs),
    MinCG = rep(NA, N.configs),
    N.DMRs = rep(NA, N.configs)
)

ix <- 1
for (pct.sig in PERCENT_SIG_RANGE){
    for (min.cg in MIN_CG_RANGE){
        print(paste0("Trying ", 100*pct.sig, "% significant and ", min.cg, " minimum CpGs"))

        dmrs <- suppressWarnings(DSS::callDMR(input.df, pct.sig = pct.sig, minCG = min.cg))

        if (is.null(dmrs)){
            summary.df[ix, ] <- c(pct.sig, min.cg, 0)
            ix <- ix + 1
            next()
        }

        print(paste0("Found ", nrow(dmrs), " DMRs"))

        dmrs.2 <- summarize_all_dmrs(dmrs)
        dmrs.3 <- add_gene_information(dmrs.2)

        dmrs.3$DMR.pct.sig <- pct.sig
        dmrs.3$DMR.minCG <- min.cg

        ofile <- make_oname(pct.sig, min.cg)
        fwrite(dmrs.3, ofile, sep="\t")

        summary.df[ix, ] <- c(pct.sig, min.cg, nrow(dmrs))
        ix <- ix + 1
    }
}

# Write out the number of DMRs for each parameter
fwrite(summary.df, file.path(ODIR, "DMRsParamSweep.csv"))

dmrs <- DSS::callDMR(input.df)
dmrs.2 <- summarize_all_dmrs(dmrs)
dmrs.3 <- add_gene_information(dmrs.2)

ofile <- file.path(ODIR, "DMRegions.bed")
fwrite(dmrs.3, ofile, sep="\t")
#END