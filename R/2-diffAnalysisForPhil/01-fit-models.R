# 01-find-DMRs.R
# Runs Dispersion Shrinkage Estimation method on methylation data (one chromosome at a time)
suppressPackageStartupMessages({
    library(data.table)
    library(tidyverse)
    library(argparse)
    library(DSS)
})

# setwd("methylome-diff-analysis/R/2-diffAnalysisForPhil/")


parser <- ArgumentParser()
# Things that may will change
parser$add_argument("--idir", default= "../../dataDerived/methylBedImputedSubset1/withoutMCI/", help='Directory to run DSS on')
parser$add_argument("--chr", default= "chr22", help='Chromosome to run DSS on')
parser$add_argument("--covariates", default= "diagnostic_group_coded,has_risk_allele,bmi,age_at_visit,Gran,CD8T,CD4T,NK,Bcell")
parser$add_argument("--test_covariate", default= "has_risk_allele")
parser$add_argument("--num_pcs", default= 2, help= 'Number of principal components to include in analysis')
parser$add_argument("--smoothing", default= 150, help= 'Width of smoothing window')
args <- parser$parse_args()

covariates <- as.vector(unlist(str_split(args$covariates, ",")))
if (args$num_pcs > 0){
    covariates <- c(covariates, paste0("PC", 1:args$num_pcs))
}
dss.formula <- formula(paste0("~", paste(covariates, collapse = "+")))

#TODO--send this code to generate master... don't want to have to repeat this
master.df <- read_csv("../../dataDerived/masterSamplesheet.csv", show_col_types=F) %>%
    mutate(has_risk_allele = ifelse(apoe_e1 == 4 | apoe_e2 == 4, 1, 0)) %>%
    mutate(diagnostic_group_coded = ifelse(diagnostic_group == "CONTROL", 0, 
            ifelse(diagnostic_group == "MCI", 1, 2)))


#Make sure the numbers are sensible of APOE risk vs protective
master.df %>%
    group_by(has_risk_allele, apoe_e1, apoe_e2) %>%
    summarize(count = n()) 

#Load data with a bit of munging. Column names are better standardized on the full dataset
blood.df <- fread(file.path(args$idir, "blood-cell-composition.csv")) %>% 
    dplyr::rename(sample_id = "sample") %>% arrange(sample_id)
PCs <- fread(file.path(args$idir, paste0(args$chr, ".PCs.csv"))) %>% arrange(sample_id)

#M(ethylated) reads and Cov(erage) from sequencing
M <- fread(file.path(args$idir, paste0(args$chr, ".M.bed")))
Cov <- fread(file.path(args$idir, paste0(args$chr, ".Cov.bed")))

df <- full_join(blood.df, master.df, by = "sample_id") %>%
        full_join(PCs, by = "sample_id") %>%
        select(c("sample_id", "diagnostic_group", all_of(covariates))) %>%
        mutate(diagnostic_group = factor(diagnostic_group, levels = c("CONTROL", "LOAD"))) %>%
        dplyr::filter(sample_id %in% names(M)[-1])

#Impute one(?) missing BMI
df$bmi[is.na(df$bmi)] <- mean(df$bmi, na.rm = TRUE)

#Check that ordering of samples in M/Cov is the same as in samplesheet
#should be in ascending order
if (all(df$sample == names(M)[-1])){
    print("Sample order is correct...")
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
test.var <- "has_risk_allele"
test.result <- DMLtest.multiFactor(dml.fit, coef = test.var)

#To measure effect size, other beta coefficients...
beta.df <- data.frame(dml.fit$fit$beta)
colnames(beta.df) <- colnames(dml.fit$X)

#Save the models
#TODO: odir based on test-covariate
odir <- file.path(args$idir, "tests", args$test_covariate)
dir.create(odir, showWarnings = F)
outname <- file.path(odir, paste0(args$chr, ".models.RData"))

save(list = intersect(ls(), c("beta.df", "test.result", "df", "dss.formula", "test.var")), file = outname)
