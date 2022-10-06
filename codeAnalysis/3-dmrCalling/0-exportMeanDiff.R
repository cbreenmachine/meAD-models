suppressPackageStartupMessages({
    library(data.table)
    library(tidyverse)
    library(argparse)
    # library(genomation)
})


parser <- ArgumentParser(description='')
parser$add_argument('--idir', default = "../../dataDerived/methylBedImputed-controlLOAD/", help = "input directory.")
parser$add_argument('--master_file', default = "../../dataDerived/masterSamplesheet.csv", help = "input directory.")
parser$add_argument('--odir', default = "../../dataSummaries/", help = "output directory.")
args <- parser$parse_args()


df <- read_csv(args$master_file, show_col_types=F)
valid_cols <- names(fread(file.path(args$idir, "chr1.M.bed"), nrow=1))
df.2 <- df %>% dplyr::filter(sample_id %in% valid_cols)


load.samples <- as.character(df.2$sample_id[df.2$diagnostic_group == "LOAD"])
control.samples <- as.character(df.2$sample_id[df.2$diagnostic_group == "CONTROL"])


get_mean_diff <- function(dir, chr, load.samples, control.samples){

    M <- fread(file.path(dir, paste0(chr, ".M.bed")))
    Cov <- fread(file.path(dir, paste0(chr, ".Cov.bed")))

    me <- M[, -"chromStart"] / Cov[ ,-"chromStart"]
    LOAD.minus.CONTROL <- rowMeans(me[, ..load.samples]) - rowMeans(me[, ..control.samples])
    chrom <- rep(chr, nrow(me))

    return( cbind(chrom, M[ ,'chromStart'], LOAD.minus.CONTROL))
}

wrapper <- function(x){
    return(get_mean_diff(args$idir, chr=x, load.samples, control.samples))
}


all_chroms <- paste0("chr", 1:22)

diff.df <- do.call("rbind", lapply(all_chroms, wrapper)) 

diff.df$chromEnd <- diff.df$chromStart + 1

ofile <- file.path(args$odir, "loadMinusControl.bed")
fwrite(diff.df, ofile)


