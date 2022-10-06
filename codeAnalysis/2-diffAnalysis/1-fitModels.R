# 01-find-DMRs.R
# Runs Dispersion Shrinkage Estimation method on methylation data (one chromosome at a time)
suppressPackageStartupMessages({
    library(data.table)
    library(tidyverse)
    library(argparse)
    library(DSS)
})

# setwd("/z/Comp/kelesgroup/loadmethyseq/methylome-diff-analysis/R/2-diffAnalysis")
parser <- ArgumentParser()
parser$add_argument("--idir", default= "../../dataDerived/methylBedImputed-controlMCILOAD/", help='Directory to run DSS on')
parser$add_argument("--odir", default= "../../dataDerived/model-controlMCILOAD/", help='Directory to run DSS on')
parser$add_argument("--config", default="controlMCILOAD")

parser$add_argument("--blood_file", default="../../dataDerived/methylBedForDXM/bloodCellPropsDXM.csv")
parser$add_argument("--master_file", default="../../dataDerived/masterSamplesheet.csv")
parser$add_argument("--chr", default= "chr22", help='Chromosome to run DSS on')
parser$add_argument("--covariates", default= "diagnostic_group_coded,bmi,sex,age_at_visit,type1,type2,type3,type4,type5")

parser$add_argument("--test_covariate", default= "diagnostic_group_coded")
parser$add_argument("--num_pcs", default= 2, help= 'Number of principal components to include in analysis')
parser$add_argument("--smoothing", default= 150, help= 'Width of smoothing window')
args <- parser$parse_args()


# config_matches_input <- str_detect(args$idir, args$config)
# config_matches_output <- str_detect(args$odir, args$config)

# if (config_matches_input & config_matches_output){
#     print("Input directory, output directory, and configuration seem to match...")
# } else {
#     stop("idir, odir don't match with configuration...")
# }

# args$include_MCI = TRUE
# args$idir = "../../dataDerived/methylBedImputed/withMCI/"
#Logic of arguments --------------------------------------------------------------------------------------------
#Cobble "baseline" covariates with PCs
covariates <- as.vector(unlist(str_split(args$covariates, ",")))
if (args$num_pcs > 0){
    covariates <- c(covariates, paste0("PC", 1:args$num_pcs))
}
dss.formula <- formula(paste0("~", paste(covariates, collapse = "+")))

#TODO--send this code to generate master... don't want to have to repeat this
master.df <- read_csv(args$master_file, show_col_types=F) %>%
    arrange(sample_id) %>% 
    dplyr::mutate(sample_id = as.character(sample_id)) %>%
    dplyr::mutate(diagnostic_group_coded = ifelse(diagnostic_group == "CONTROL", 0, 
        ifelse(diagnostic_group == "MCI", 0.5, 1)))


#Load data with a bit of munging. Column names are better standardized on the full dataset
blood.df <- read_csv(args$blood_file, show_col_types=F) %>% 
    arrange(sample_id) %>% 
    dplyr::mutate(sample_id = as.character(sample_id))

PCs <- read_csv(file.path(args$idir, paste0(args$chr, ".PCs.csv")), show_col_types=F) %>% 
    arrange(sample_id) %>% 
    dplyr::mutate(sample_id = as.character(sample_id))

#M(ethylated) reads and Cov(erage) from sequencing
M <- fread(file.path(args$idir, paste0(args$chr, ".M.bed")))
Cov <- fread(file.path(args$idir, paste0(args$chr, ".Cov.bed")))

# valid_sample_ids <- intersect(names(M), df$sample_id)

#Pull phenotypes, PCs, and blood into one design matrix
#Subset to only have sample_ids also contained in M/Cov
df <- full_join(blood.df, master.df, by = "sample_id") %>%
        full_join(PCs, by = "sample_id") %>%
        select(c("sample_id", all_of(covariates))) %>%
        drop_na() %>%
        distinct()

if (args$config == "controlMCILOAD"){df <- dplyr::filter(df, diagnostic_group_coded != 1/2)}

# Recompute in case any NAs caused us to drop a sample
valid_sample_ids <- intersect(names(M), df$sample_id)

df <- dplyr::filter(df, sample_id %in% valid_sample_ids)

keep_cols <- c("chromStart", sort(as.character(valid_sample_ids)))
M <- M %>% dplyr::select(all_of(keep_cols))
Cov <- Cov %>% dplyr::select(all_of(keep_cols))

print(paste0("Number of samples : ", nrow(df)))

#Check that ordering of samples in M/Cov is the same as in samplesheet
#should be in ascending order
if (all(df$sample_id == names(M)[-1])){
    print("Sample order is correct...")
} else {
    warning("Sample order wrong!")
}

# create bs seq object, needs chromosome identiifer, methylated reads, and unmethylated reads
bs <- BSseq(chr = rep(args$chr, nrow(M)), pos = M$chromStart,
            M = as.matrix(M[ , -c("chromStart"), with=FALSE]), 
            Cov = as.matrix(Cov[, -c("chromStart"), with=FALSE]), 
            sampleNames = names(M)[-1])

# Derive some parameters
smooth = TRUE
if (args$smoothing == 0){
    smooth = FALSE
}


dml.fit <- DMLfit.multiFactor(bs, design = df, smoothing = smooth, 
    smoothing.span = args$smoothing, formula = dss.formula)


#TODO: interpret this from command line args...
test.var <- args$test_covariate
test.result <- DMLtest.multiFactor(dml.fit, coef = test.var)

#To measure effect size, other beta coefficients...
beta.df <- data.frame(dml.fit$fit$beta)
colnames(beta.df) <- colnames(dml.fit$X)

#Save the models
#TODO: odir based on test-covariate
odir <- file.path(args$odir, args$test_covariate)
dir.create(odir, showWarn=F, recursive=T)
outname <- file.path(odir, paste0(args$chr, ".models.RData"))

save(list = intersect(ls(), c("beta.df", "test.result", "df", "dss.formula", "test.var")), file = outname)
