# 02-call-DMRs.R
# Runs Dispersion Shrinkage Estimation method on methylation data (one chromosome at a time)
suppressPackageStartupMessages({
    library(data.table)
    library(tidyverse)
    library(fdrtool)
    library(argparse)
    library(DSS)
})


N <- 84
n_cov <- 12
degrees_freedom <- N - n_cov + 1

parser <- ArgumentParser()
parser$add_argument("--dir", default= "../../data/07-counts-by-chrom-imputed-subset2/", help='Directory to run DSS on')
args <- parser$parse_args()

all.files <-list.files(file.path(args$dir, "models"), full=T)


# test.cohort has chr,pos,stat,pvals,fdrs as columns and is DSS-sepcific instance
.load_data <- function(path){
    # Load from path like "./chr8/smooth-150-PCs-2/models.RData"
    load(path)
    df <- cbind(test.cohort, beta.df) %>% drop_na() 

    # This allows us to do testing with DSS framework
    class(df)[2] <- "DMLtest.multiFactor"
    return(df)
}


.split_by_chrom <- function(df){
    out <- list()

    for (ii in 1:22){
        out[[ii]] <- df[df$chr == paste0("chr", ii), ]
    }
    return(out)
}


.correct_pvals <- function(df, deg=degrees_freedom){
    # Leaves other columns in tact
    df$DSS.stat <- df$stat
    df$DSS.pvals <- df$pvals
    df$DSS.fdrs <- df$fdrs

    # Extract test statistic and adjust
    zz <- (df$DSS.stat - mean(df$DSS.stat)) / sd(df$DSS.stat)
    out <- fdrtool(zz, statistic = "normal", plot = FALSE)

    # Overrides raw values from DSS
    df$pvals <- out$pval
    df$lfdrs <- out$lfdr 
    df$stat <- qchisq(out$pval, df=deg, lower.tail=FALSE)

    # New values
    df$qvals <- out$qval
    
    return(df)
}


# How to do this efficiently... Loop thru chromosomes??
.load_wrapper <- function(path){ .load_data(path) %>% return()}
df <- do.call("rbind", lapply(all.files, .load_wrapper)) 

.correct_wrapper <- function(df){ .correct_pvals(df) %>% return()}
df.2 <- do.call("rbind", lapply(.split_by_chrom(df), .correct_wrapper)) 

# Call DMRs on both...
dmrs.raw <- DSS::callDMR(df)
dmrs.corr  <- DSS::callDMR(df.2)

dmrs.raw %>% transmute(
    chrom = chr, 
    chromStart = start - 1, 
    chromEnd = end + 1,
    length, nCG, areaStat
    ) %>%
    fwrite(file.path(args$dir, "DMRs.raw.bed"), sep="\t")


dmrs.corr %>% transmute(
    chrom = chr, 
    chromStart = start - 1, 
    chromEnd = end + 1,
    length, nCG, areaStat
    ) %>%
    fwrite(file.path(args$dir, "DMRs.adj.bed"), sep="\t")


class(df.2) <- "data.frame"
DMPs.df <- df.2 %>% 
    dplyr::mutate(chrom = chr, chromStart = pos - 1, chromEnd = chromStart + 1) %>% 
    dplyr::select(-pos) 

ofile <- file.path(args$dir, "DMPs.bed")
fwrite(DMPs.df, ofile, sep="\t")

#END