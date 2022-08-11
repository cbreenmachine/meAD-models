#/ua/cebreen/miniconda3/envs/R/bin/R
# Runs Dispersion Shrinkage Estimation method on methylation data (one chromosome at a time)
suppressPackageStartupMessages({
    library(data.table)
    library(tidyverse)
    library(argparse)
    library(methylCC)
    library(GenomicRanges)
    library(minfi)
    library(bsseq)
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
parser$add_argument("--chr", default= "chr1", help='Chromosome to impute on')
parser$add_argument("--idir", default= "../../data/07-counts-by-chrom-imputed-subset1/")
args <- parser$parse_args()

.read_data <- function(s=".M.bed"){
    out <- fread(file.path(args$idir, paste0(args$chr, s))) %>%
            dplyr::rename(start = chromStart) %>%
            dplyr::mutate(chrom = args$chr, end = start + 1) %>%
            GRanges()
    return(out)
}

M.gr <- .read_data()
Cov.gr <- .read_data(".Cov.bed")

# N <- nrow(M)

# M.gr <- GenomicRanges::GRanges(dplyr::rename(data.frame(M), start = chromStart))

# bs <- BSseq(M = as.matrix(M[, -"chromStart"]), Cov = as.matrix(Cov[, -"chromStart"]), 
#             pos = as.vector(M$chromStart) + 1, chr = rep(args$chr, N))

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

bs <- lift_to_hg19(M.gr, Cov.gr, chain)

est <- estimatecc(object = bs, include_cpgs = TRUE, include_dmrs=TRUE) 
