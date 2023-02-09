# 01-find-DMRs.R
# Runs Dispersion Shrinkage Estimation method on methylation data (one chromosome at a time)
suppressPackageStartupMessages({
    library(argparse)
    library(DSS)
    library(magrittr)
    library(dplyr)
    library(stringr)
    library(parallel)
})

parser <- ArgumentParser()
parser$add_argument("--idir", default= "../DataDerived/ControlLOAD/Inputs/")
parser$add_argument("--odir", default= "../DataDerived/ControlLOAD/test-diagnostic-group-coded/")
parser$add_argument("--chr", default= "chr22")
parser$add_argument("--baseline_covariates", default= "bmi,sex,age_at_visit,type1,type2,type3,type4,type5,PC1,PC2")
parser$add_argument("--test_covariate", default= "diagnostic_group_coded")
parser$add_argument("--smoothing", default= 150, help= 'Width of smoothing window')
args <- parser$parse_args()

ifile <- file.path(args$idir, paste0("filtered-DSS-inputs-", args$chr, ".RData"))

###########################################################
##################### Directory Manipulation ##############
###########################################################
load(ifile)

# Temp patch, fixed in 0-prepare
# Should not have been storing a `args`
args <- parser$parse_args()

# Put the two paths together and create the output
odir <- args$odir
dir.create(odir, showWarn=F, recursive=T)

file_name <- str_remove(str_replace(basename(ifile), "input", "output"), "filtered-")
outname <- file.path(odir, file_name)

###########################################################
##################### DATA PREP ###########################
###########################################################
covariates_split <- stringr::str_split(args$baseline_covariates, ",")[[1]]
dss.formula <- formula(paste0("~", args$test_covariate, "+", paste(covariates_split, collapse = "+")))

# Pull out just the covariates we'll use
design.df <- df %>% 
    tibble::column_to_rownames("sample_id") %>%
    dplyr::select(all_of(c(args$test_covariate, covariates_split))) %>%
    dplyr::mutate_at(args$test_covariate, as.numeric) %>% 
    tidyr::drop_na()

# If there are NAs in any of the covariates specified, be sure to drop
#corresponding sample in bs object
ix <- colnames(bs) %in% rownames(design.df)
bs.sub <- bs[ , ix]

if (all(colnames(bs.sub) == rownames(design.df))){
    print("After filtering NAs, rownames in design match colnames in bs")
} else {
    warning("rownames in design != colnames in bs")
}

# Derive some parameters
smooth = TRUE
if (args$smoothing == 0){
    smooth = FALSE
}

# Fit models
dml.fit <- DMLfit.multiFactor(bs.sub, design = design.df, smoothing = smooth, 
    smoothing.span = args$smoothing, formula = dss.formula)

invisible(gc())

# Test covariate (to start this will be LOAD status--diagnostic_group_coded)
test.var <- args$test_covariate
test.result <- DMLtest.multiFactor(dml.fit, coef = test.var)

#To measure effect size, other beta coefficients...
beta.df <- data.frame(dml.fit$fit$beta)
colnames(beta.df) <- colnames(dml.fit$X)

design.mat <- model.matrix(dss.formula, design.df)

save(list = intersect(ls(), c("beta.df", "test.result", "design.df", "dss.formula", "test.var", "design.df")), 
    file = outname)

#END