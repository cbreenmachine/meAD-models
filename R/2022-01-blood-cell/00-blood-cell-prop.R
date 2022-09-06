#/ua/cebreen/miniconda3/envs/R/bin/R

suppressPackageStartupMessages({
    library(data.table)
    library(tidyverse)
    library(argparse)
    library(methylCC)
    library(GenomicRanges)
    library(minfi)
    library(bsseq)
    library(FlowSorted.Blood.450k)

})


chain.file <- "hg38ToHg19.over.chain"
.check_chain <- function(){
    if (! file.exists(chain.file)){
        tmp <- file.path("https://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/", chain.file)
        utils::download.file(tmp, destfile = paste0(chain.file, ".gz"))
        system(paste0("gzip -d ", chain.file, ".gz"))
    } 
}
.check_chain()
chain <- rtracklayer::import.chain(chain.file)


parser <- ArgumentParser()
parser$add_argument("--idir", default= "../../data/07-counts-by-chrom-imputed-subset1/")
args <- parser$parse_args()

# BED files are zero-based 
read_data <- function(idir, chrom, s) {
    out <- fread(file.path(idir, paste0(chrom, ".", s, ".bed"))) %>%
            dplyr::rename(start = chromStart) %>%
            dplyr::mutate(chrom = chrom, start = start) %>%  #TODO: bed file is not correctly positioned
            dplyr::mutate(end = start + 1) %>%
            makeGRangesFromDataFrame(starts.in.df.are.0based = T, keep.extra.columns = T)
    return(out)
}



# Helper function for lifting...
clean <- function(x){return(str_remove(x, "X"))}


lift_to_hg19 <- function(M.gr, Cov.gr, chain){
    M.lifted <- rtracklayer::liftOver(M.gr, chain) %>% data.frame()
    Cov.lifted <- rtracklayer::liftOver(Cov.gr, chain)
    
    M  <- data.frame(M.gr) %>% rename_with(clean)
    Cov <- data.frame(Cov.gr) %>% rename_with(clean)

    bs <- BSseq(M = as.matrix(M[, -1:-5]), Cov = as.matrix(Cov[, -1:-5]), 
        pos = as.vector(M$start), chr = M$seqnames)
    return(bs)
}


read_wrapper <- function(chr){
    Cov <- read_data(args$idir, chr, "Cov")
    M <- read_data(args$idir, chr, "M")
    bs <- lift_to_hg19(M, Cov, chain)
    return(bs)
}


all_chromosomes <- paste0("chr", seq(1,22))
all_bs_objects <- lapply(all_chromosomes, read_wrapper)

bs <- do.call(rbind, all_bs_objects)

est <- estimatecc(object = bs, include_cpgs=TRUE, include_dmrs=TRUE)


df_methylCC = gather(cbind("samples" = rownames(cell_counts(est)),
                           cell_counts(est)),
                     celltype, est, -samples)