# 01-find-DMRs.R
# Runs Dispersion Shrinkage Estimation method on methylation data (one chromosome at a time)
suppressPackageStartupMessages({
    library(data.table)
    library(tidyverse)
    library(argparse)
    library(DSS)
})

parser <- ArgumentParser()
# Things that may will change
parser$add_argument("--chr", default= "chr22", help='Chromosome to run DSS on')
parser$add_argument("--covariates", default= "Gran,CD8T,CD4T,NK,Bcell,bmi,age,sex")
parser$add_argument("--num_pcs", default= 2, help= 'Number of principal components to include in analysis')
parser$add_argument("--smoothing", default= 150, help= 'Width of smoothing window')
# File paths
parser$add_argument("--samples_file", default= "../../data/study-data/array-samples.csv", help="CSV file with samples to filter")
parser$add_argument("--pheno_file", default= "../../data/study-data/phenos-cleaned.csv", help="CSV file with LOAD/control status")
parser$add_argument("--pc_dir", default= "../../data/old-data/prin-comps-array-samples/", help="CSV file with LOAD/control status")
args <- parser$parse_args()

covariates <- as.vector(unlist(str_split(args$covariates, ",")))
if (args$num_pcs > 0){
    covariates <- c(covariates, paste0("PC", 1:args$num_pcs))
}
dss.formula <- as.formula(paste(c("~cohort", covariates), collapse = "+"))


idir <- file.path("./old-data/", args$chr)
odir <- paste0(idir, "/smooth-", as.character(args$smoothing), "-PCs-", args$num_pcs, "/")
dir.create(odir, showWarnings = FALSE)

pc.df <- read_csv(file.path(args$pc_dir,  paste0(args$chr, ".csv")), col_types = cols())
ss.df <- read_csv(args$pheno_file, col_types = cols())

M <- fread(file.path(idir, "M.csv"))
Cov <- fread(file.path(idir, "Cov.csv"))

#TODO: restructure 71 to look like subset
#TODO: just take one directory as input, and chrom

valid.samples <- intersect(intersect(colnames(M), ss.df$sample), pc.df$sample)

filt.df <- ss.df %>%
            inner_join(pc.df, by = "sample") %>% 
            dplyr::filter(sample %in% valid.samples) 

drop.cols <- as.character(setdiff(colnames(M), valid.samples))[-1] # dont remove pos

if (length(drop.cols) > 0){
    M[, (drop.cols):=NULL]
    Cov[, (drop.cols):=NULL]
}

#--> Handle missing values?
#filt.df$pack_years[is.na(filt.df$pack_years)] <- 0

# Handle mismatching files
all(filt.df$sample == names(M)[-1])

# create bs seq object, needs chromosome identiifer, methylated reads, and unmethylated reads
bs <- BSseq(chr = rep(args$chr, nrow(M)), pos = M$pos,
            M = as.matrix(M[ , -c("pos"), with=FALSE]), 
            Cov = as.matrix(Cov[, -c("pos"), with=FALSE]), 
            sampleNames = names(M)[-1])

all( filt.df$sample == colnames(bs) )


# Derive some parameters
smooth = TRUE
if (args$smoothing == 0){
    smooth = FALSE
}


dml.fit <- DMLfit.multiFactor(bs, design = filt.df, smoothing = smooth, 
smoothing.span = args$smoothing, formula = dss.formula)

test.var <- colnames(dml.fit$X)[2]
print(test.var)

print("First five in targets: ")
filt.df$cohort[1:5]

print("First five in cohortCONTROL:")
as.vector(dml.fit$X[1:5, test.var])

# In cohortCONTROL....
# CONTROL coded as 1
# LOAD coded as 0
# beta represents going "from LOAD to CONTROL"
# positive beta means control methylation is higher (LOAD is lower) (so hypomethylated)
# Need to flip sign to "make sense"

test.cohort <- DMLtest.multiFactor(dml.fit, coef = test.var)

#--> To measure effect size, other beta coefficients...
beta.df <- data.frame(dml.fit$fit$beta)
colnames(beta.df) <- colnames(dml.fit$X)
beta.df <- beta.df %>% mutate(LOAD.minus.control = - cohortCONTROL)

# Save the models
outname <- file.path(odir, "models.RData")
print(outname)
save(list = intersect(ls(), c("beta.df", "test.cohort", "filt.df")), file = outname)
