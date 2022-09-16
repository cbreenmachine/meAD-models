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
parser$add_argument("--dir", default= "../../data/07-counts-by-chrom-imputed-subset1/", 
                    help='Directory to run DSS on')
parser$add_argument("--chr", default= "chr22", help='Chromosome to run DSS on')
parser$add_argument("--covariates", default= "Gran,CD8T,CD4T,NK,Bcell,BMI,age_at_visit,sex")
parser$add_argument("--num_pcs", default= 2, help= 'Number of principal components to include in analysis')
parser$add_argument("--smoothing", default= 150, help= 'Width of smoothing window')
args <- parser$parse_args()

covariates <- as.vector(unlist(str_split(args$covariates, ",")))
if (args$num_pcs > 0){
    covariates <- c(covariates, paste0("PC", 1:args$num_pcs))
}
dss.formula <- as.formula(paste(c("~diagnostic_group", covariates), collapse = "+"))


#--> Data loading
blood.df <- fread(file.path(args$dir, "blood-cell-composition.csv"))
master.df <- fread(file.path(args$dir, "master.csv"))
PCs <- fread(file.path(args$dir, paste0(args$chr, ".PCs.bed")))


df <- full_join(blood.df, master.df, by = c("sample" = "RA_id")) %>%
        full_join(PCs, by = c("sample" = "sample_id")) %>%
        select(c("sample", "diagnostic_group", all_of(covariates))) %>%
        mutate(diagnostic_group = factor(diagnostic_group, levels = c("Control", "LOAD"))) %>%
        column_to_rownames("sample")


#--> Missing BMI...
# cell_types <- c("Gran", "CD8T", "CD4T", "NK", "Bcell")
df$BMI[is.na(df$BMI)] <- mean(df$BMI, na.rm = TRUE)

M <- fread(file.path(args$dir, paste0(args$chr, ".M.bed")))
Cov <- fread(file.path(args$dir, paste0(args$chr, ".Cov.bed")))

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

test.var <- colnames(dml.fit$X)[2]
print(test.var)

print("First five in targets: ")
df$diagnostic_group[1:5]

print("First five in cohortCONTROL:")
as.vector(dml.fit$X[1:5, test.var])

# In cohortCONTROL....
# CONTROL coded as 0
# LOAD coded as 1
# beta represents going "from CONTROL to LOAD"
test.cohort <- DMLtest.multiFactor(dml.fit, coef = test.var)

#--> To measure effect size, other beta coefficients...
beta.df <- data.frame(dml.fit$fit$beta)
colnames(beta.df) <- colnames(dml.fit$X)
beta.df <- beta.df %>% mutate(LOAD.minus.control = diagnostic_groupLOAD)

# Save the models
odir <- file.path(args$dir, "models")
dir.create(odir, showWarnings = F)

outname <- file.path(odir,  paste0(args$chr, "-models.RData"))
print(outname)

save(list = intersect(ls(), c("beta.df", "test.cohort", "df")), file = outname)
