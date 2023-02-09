
# Read in one by one the pi files and subset by overlaps
# 5-subset-pi-by-dmr.R
suppressPackageStartupMessages({
    library(data.table)
    library(magrittr)
    library(argparse)
    library(GenomicRanges)
    library(bsseq)
})


parser <- ArgumentParser()
parser$add_argument(
    "--idir", 
    default= "../DataDerived/ControlLOAD/Inputs/", 
    help='Where the input-chrX.RData files are stored')

parser$add_argument(
    "--dmr_file", 
    default= "../DataDerived/ControlLOAD/test-diagnostic-group-coded/Summaries/DMRegions.bed", 
    help='path to DMR file')
args <- parser$parse_args()

# Read DMRs and cast to GRanges
dmrs <- fread(args$dmr_file) %>%
    makeGRangesFromDataFrame(starts = T, keep.extra = T)
end(dmrs) <- end(dmrs) + 1 # DSS doesn't know


cast_to_granges <- function(DT){
    makeGRangesFromDataFrame(DT, starts = T, keep.extra = T)
}

subset_bs_by_ix <- function(bs, ix, tp){
    # tp is "M" or "Cov"

    out <- as.data.frame(getCoverage(bs[ix, ], type = tp))
    out$chr <- as.vector(seqnames(bs[ix, ]))
    out$start <- as.vector(start(bs[ix, ]))
    out$end <- as.vector(end(bs[ix, ]))
    out$type <- tp
    dplyr::select(out, c("chr", "start", "end", "type", everything()))
}

pipeline <- function(ff, regions=dmrs){
    load(ff)
    print(ff)

    # Make it 
    end(bs) <- start(bs) + 2

    # medians <- apply(getCoverage(bs), 1, median)

    # Find overlaps, extract coordinates
    overlaps <- findOverlaps(bs, dmrs, min=1)

    posix <- queryHits(overlaps)
    regix <- subjectHits(overlaps)

    # Subset region (duplicate rows)
    sub <- as.data.frame(dmrs[regix,])

    # Rename to avoid duplicate chr, start, end
    names(sub)[1:3] <- c("chr.dmr", "start.dmr", "end.dmr")

    M.df <- cbind(sub, subset_bs_by_ix(bs, posix, "M"))
    Cov.df <- cbind(sub, subset_bs_by_ix(bs, posix, "Cov"))

    # Return as one object
    rbind(M.df, Cov.df)
}


# all.files <- list.files(args$idir, pattern = "input*", full=T)
tmp <- as.vector(unique(seqnames(dmrs)))
all.files <- file.path(args$idir, paste0("filtered-DSS-inputs-", tmp, ".RData"))

out <- do.call(rbind, lapply(all.files, FUN=pipeline))
ofile <- file.path(dirname(args$dmr_file), "DMRegions-CpGs.bed")
fwrite(out, ofile, sep = "\t")

#END