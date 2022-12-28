
# Read in one by one the pi files and subset by overlaps
# 5-subset-pi-by-dmr.R
suppressPackageStartupMessages({
    library(data.table)
    library(magrittr)
    library(argparse)
    library(GenomicRanges)
})


parser <- ArgumentParser()
parser$add_argument(
    "--idir", 
    default= "../DataDerived/ControlLOAD/test-diagnostic-group-coded/", 
    help='Where the pi files are stored')

parser$add_argument(
    "--dmr_file", 
    default= "../DataDerived/ControlLOAD/test-diagnostic-group-coded/ExperimentSummary/DMRegions-filtered.bed", 
    help='path to DMR file')
args <- parser$parse_args()

# Read DMRs and cast to GRanges
dmrs <- fread(args$dmr_file) %>%
    makeGRangesFromDataFrame(starts = T, keep.extra = T)
end(dmrs) <- end(dmrs) + 1 # DSS doesn't know


all.files <- list.files(args$idir, pattern = "pi\\.chr*", full=T)

cast_to_granges <- function(DT){
    makeGRangesFromDataFrame(DT, starts = T, keep.extra = T)
}

pipeline <- function(ff, regions=dmrs){
    DT <- fread(ff)
    GR <- cast_to_granges(DT)

    # Find overlaps, extract coordinates
    overlaps <- findOverlaps(GR, dmrs, min = 1)

    posix <- queryHits(overlaps)
    regix <- subjectHits(overlaps)

    # Subset region (duplicate rows)
    sub <- as.data.frame(dmrs[regix,])

    # Rename to avoid duplicate chr, start, end
    names(sub)[1:3] <- c("chr.dmr", "start.dmr", "end.dmr")

    # Return as one object
    cbind(sub, DT[posix, ])
}

out <- do.call(rbind, lapply(all.files, FUN=pipeline))

ofile <- file.path(dirname(args$dmr_file), "DMR.pis.bed")
fwrite(out, ofile)

# What's going on at these positions?
