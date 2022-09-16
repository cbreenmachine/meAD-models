# 02-call-DMRs.R
# Runs Dispersion Shrinkage Estimation method on methylation data (one chromosome at a time)
suppressPackageStartupMessages({
    library(data.table)
    library(tidyverse)
    library(fdrtool)
    library(argparse)
    library(DSS)
})

parser <- ArgumentParser()
parser$add_argument("--idir", default= "../../dataDerived/methylBedImputed/withoutMCI/tests/diagnostic_group_coded/", help='Where are the models stored')
parser$add_argument("--odir", default= "../../dataDerived/DMPs/noMCI-withBlood", help='Where are the models stored')
args <- parser$parse_args()

dir.create(args$odir, recursive=T, showWarn=F)
all.files <-list.files(args$idir, pattern = "*RData", full=T)

#Derive the degrees of freedom for p-value correction
load(all.files[1])
N <-nrow(df)
n_cov <- length(names(beta.df))
degrees_freedom <- N - n_cov + 1

# test.cohort has chr,pos,stat,pvals,fdrs as columns and is DSS-sepcific
.load_data <- function(path){
    # Load from path like "./chr8/smooth-150-PCs-2/models.RData"
    load(path)
    df <- cbind(test.result, beta.df) %>% drop_na() 

    # This allows us to do testing with DSS framework
    class(df)[2] <- "DMLtest.multiFactor"
    return(df)
}


.split_by_chrom <- function(df){
    out <- list()

    valid_nums <- as.numeric(str_replace(unique(df$chr), "chr", ""))
    ix <- 1

    for (ii in valid_nums){
        out[[ix]] <- df[df$chr == paste0("chr", ii), ]
        ix <- ix + 1
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
# dmrs.raw <- DSS::callDMR(df)
dmrs.corr  <- DSS::callDMR(df.2)
# dmrs.raw %>% transmute(
#     chrom = chr, 
#     chromStart = start - 1, 
#     chromEnd = end + 1,
#     length, nCG, areaStat
#     ) %>%
#     fwrite(file.path(args$dir, "DMRs.raw.bed"), sep="\t")


dmrs.corr %>% transmute(
    chrom = chr, 
    chromStart = start - 1, 
    chromEnd = end + 1,
    length, nCG, areaStat
    ) %>%
    fwrite(file.path(args$dir, "DMRs.adj.bed"), sep="\t")

# CpG starts on position 5 (how it's stored right now)
# 1-based closed would say it goes from [5, 6]
# 0-based open would say it goes from (4, 6]
class(df.2) <- "data.frame"
DMPs.df <- df.2 %>% 
    dplyr::mutate(chrom = chr, chromStart = pos - 1, chromEnd = pos + 1) %>% 
    dplyr::select(-pos) %>%
    dplyr::relocate(c(chrom, chromStart, chromEnd))

ofile <- file.path(args$odir, "DMPs.bed")
fwrite(DMPs.df, ofile, sep="\t")

#END