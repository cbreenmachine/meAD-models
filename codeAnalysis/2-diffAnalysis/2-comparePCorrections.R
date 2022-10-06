# Bacon adjustment
# FDR Tool adjustment
# Lambda adjustment
# comb-p adjustment
suppressPackageStartupMessages({
    library(data.table)
    library(tidyverse)
    library(fdrtool)
    library(argparse)
    library(DSS)
    library(bacon)
})

parser <- ArgumentParser()
parser$add_argument("--idir", default= "../../dataDerived/models-controlLOAD-pc2-ct5/diagnostic_group_coded/", help='Where are the models stored')
parser$add_argument("--combp_file", default= "../../dataSummaries/dmps-controlLOAD-pc2-ct5/pvals.acf.bed", help='Where are the models stored')
parser$add_argument("--odir", default= "../../dataSummaries/pComp-controlLOAD-pc2-ct5/", help='Where are the models stored')
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
    out <- fdrtool(zz, statistic = "normal", cutoff.method = "locfdr", plot = FALSE)

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
p.dss <- df$pvals

# FDR toolkit
zz <- df$stat
out <- fdrtool(zz, statistic = "normal", cutoff.method = "locfdr", plot = FALSE)
p.fdr <- out$pval

# Bacon
bc <- bacon(df$stat, niter=20000L)
p.bacon <- pval(bc)

# Combp
combp.df <- fread(args$combp_file, sep="\t") %>% dplyr::filter(V5 > 0)
p.combp <- combp.df$V5

p_to_lambda <- function(pvals){
    chisq <- qchisq(1-pvals, 1)
    return(median(chisq) / qchisq(0.5, 1))
}


# 
lambda.dss <- p_to_lambda(p.dss)
lambda.bacon <- p_to_lambda(p.bacon)
lambda.fdr <- p_to_lambda(p.fdr)
lambda.combp <- p_to_lambda(p.combp)

out.df <- data.frame(
    scheme = c("DSS", "Bacon", "FDRtool", "CombP"),
    lambda = c(lambda.dss, lambda.bacon, lambda.fdr, lambda.combp)
)


write_csv(out.df, "../../dataSummaries/2022-10-04-lambdaComparisonsControlLOAD.csv")
