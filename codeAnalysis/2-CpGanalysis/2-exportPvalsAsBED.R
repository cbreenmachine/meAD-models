# 3-exportDMPsAsBED.R
# 
suppressPackageStartupMessages({
    library(data.table)
    library(tidyverse)
    library(argparse)
    library(DSS)
    library(fdrtool)
    library(GenomicRanges)
    library(Rcpp)
})


Rcpp::sourceCpp("./.matrix_multiply.cpp")

parser <- ArgumentParser()
parser$add_argument("--idir", default= "../../dataDerived/analysis-controlLOAD/test-diagnostic-group-coded/", help='Where are the models stored')
args <- parser$parse_args()

# Output directories
last.dir <- "experimentSummary"
odir <- file.path(args$idir, last.dir)
print(paste0("Will write pvals.df, dmrs.df, and pi.df to ", odir))
dir.create(odir, recursive=T, showWarn=F)

# Load data
all.files <-list.files(args$idir, pattern = "*RData", full=T)


compute_pi_for_cpg_i <- function(i, X, beta.mat){
    # helper function that calls KMP's matrix multiply script
    #i is CpG index
    b <- beta.mat[i, ]
    tmp <- (sin(multiply(as(X,'sparseMatrix'), as(b, 'sparseMatrix'))) + 1) / 2
    as.numeric(tmp)
}


compute_pis <- function(X, beta.mat, test.result){

    N <- nrow(beta.mat)

    out <- mclapply(FUN = function(i) {
        compute_pi_for_cpg_i(i, X, beta.mat)
        }, 
        1:N,mc.cores=12
    )

    pi.mat <- do.call(rbind, out)

    # Assign column names based on sample id
    colnames(pi.mat) <- rownames(X)

    # Add chr, start, end as leading columns
    class(test.result) <- "data.frame" # remove DSS class
    bed.cols <- dplyr::transmute(test.result, chr, start = pos, end = pos + 2)

    pi.df <- data.table(cbind(bed.cols, pi.mat), rownames = NULL)
    pi.df
}

compute_pi_diffs <- function(X, design.df, pi.df){
    
    load.samples <- rownames(X)[design.df$diagnostic_group_coded == 1]
    control.samples <- rownames(X)[design.df$diagnostic_group_coded == 0]

    pi.diff <- rowMeans(pi.df[ , ..load.samples]) - 
                rowMeans(pi.df[ , ..control.samples])

}

make_pi_oname <- function(ifile){
    root <- file.path(dirname(ifile), last.dir)  
    z <- str_replace(
            str_replace(basename(ifile), "output", "pi"), 
        "RData", "bed")
    file.path(root, z)
}


run_pi_routine <- function(ifile){
    ofile <- make_pi_oname(ifile)

    if (!file.exists(ofile)){
        load(ifile)

        # While data is loaded, compute pis and 
        X <- model.matrix(dss.formula, design.df)
        beta.mat <- as(beta.df, "matrix")

        pi.df <- compute_pis(X, beta.mat, test.result)
        print(paste0("Writing out ", ofile))

        fwrite(pi.df, ofile, sep="\t")
        gc()

    } else {
        print(paste0(ofile, " already exists, skipping"))
    }
}

lapply(FUN=run_pi_routine, all.files)

# test.cohort has chr,pos,stat,pvals,fdrs as columns and is DSS-sepcific
load_data_with_pi_diff <- function(ifile){
    # Load from path like "./chr8/smooth-150-PCs-2/models.RData"
    load(ifile)

    # Very inefficient to load back in but oh well
    pi.df <- fread(make_pi_oname(ifile))

    X <- model.matrix(dss.formula, design.df)
    pi.diffs <- compute_pi_diffs(X, design.df, pi.df)

    # covariates and p-values
    df <- cbind(test.result, beta.df) 
    df$pi.diff <- pi.diffs
   
    gc()

    return(df)
}

# How to do this efficiently... Loop thru chromosomes??
df <- do.call("rbind", lapply(all.files, load_data_with_pi_diff)) 

# Genomic inflation (test stats from DSS are N(0,1))
# https://github.com/haowulab/DSS/blob/55cb1a05738a2ed02717af50b7b52828bc6b508d/R/DML.multiFactor.R#L192
chisq <- df$stat^2
lambda <- median(chisq) / qchisq(0.5, 1)
stat.crct <- chisq / lambda

# Correct by adjusting with genomic inflation 
p.crct <- pchisq(stat.crct, df=1, lower.tail=FALSE)

# Couple hundred p-values are stored as zero
ix.0 <- which(p.crct == 0)
p.crct[ix.0] <-  min(p.crct[p.crct > 0])

# Copy original df for use with DSS DMR calling
# df.2 <- df

# # this preserves sign for DSS DMR calling
# df.2$stat <- df$stat / sqrt(lambda) 
# df.2$pvals <- p.crct
# class(df.2)[2] <- "DMLtest.multiFactor"

# # Call DMRs
# p.threshold <- 1e-5 # default, could make this smaller
# pct.sig <- 0.95 # default is 0.5, so this is conservative

# dmrs <- DSS::callDMR(df.2, p.threshold = p.threshold, pct.sig = pct.sig)

# dmrs.out <- dmrs %>%
#     dplyr::transmute(chrom = .pad_chrom(chr), chromStart = start, chromEnd = end+2,
#         length, nCG, areaStat)

# ofile <- file.path(args$odir, "dmrs.DSS.bed")

# # Write out dmrs, first few lines look like
# #DMRs called with DSS w/ default parameters unless otherwise specified
# #P value threshold: 1e-05
# #Percent significant: 0.95
# # chrom	chromStart	chromEnd	length	nCG	areaStat
# # chr02	159294967	159295097	131	10	-147.077920660184
# fwrite(dmrs.out, ofile, sep="\t", append = T, col.names = T)


##############################################################
################ FDR CACLULATION #############################
##############################################################
# False discovery rate, and 
fdr.out <- fdrtool(p.crct, statistic = "pvalue", plot = FALSE)

# Rename columns, but don't sort yet
pvals.df <- df %>% 
    dplyr::mutate(chr = chr, start = pos , end = pos + 2) %>% 
    dplyr::select("chr", "start", "end", 7, "pi.diff", "stat")

# Add in the pval-type stuff we want
pvals.df$p.corrected <- p.crct
pvals.df$lfdr <- fdr.out$lfdr
pvals.df$p.raw <- df$pvals

ofile <- file.path(odir, "pvals.bed")
fwrite(pvals.df, ofile, sep="\t")
