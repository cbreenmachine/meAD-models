# 0_prepare_DSS_chrom.R
# takes M and Cov matrices, alongside phenotypes and outputs an R data store with bs object (input for DSS)
# and design matrix (df)

suppressPackageStartupMessages({
    library(data.table)
    library(argparse)
    library(DSS)
    library(magrittr)
})

parser <- ArgumentParser()
parser$add_argument("--idir", default= "../../dataDerived/methylBedImputed-controlMCILOAD/", help='Directory to run DSS on')
parser$add_argument("--blood_file", default="../../dataDerived/methylBedForDXM/bloodCellPropsDXM.csv")
parser$add_argument("--master_file", default="../../dataDerived/masterSamplesheet.csv")
parser$add_argument("--chr", default= "chr22", help='Chromosome to run DSS on')
parser$add_argument("--coverage_threshold", default=5, help="Keep sites with at least () median coverage")

args <- parser$parse_args()

read_and_arrange <- function(s){
    # s should be "M" or "Cov"
    file_name <- file.path(args$idir, paste0(args$chr, ".", s, ".bed"))
    data <- fread(file_name) %>% dplyr::arrange(chromStart)
    data
}

# Derived the output directory from the name of the input
odir <- stringr::str_replace(args$idir, "methylBedImputed", "analysis")
dir.create(odir, showWarn = F, recursive = T)
ofile <- file.path(odir, paste0("input.", args$chr, ".RData"))


#TODO--send this code to generate master... don't want to have to repeat this
master.df <- read.csv(args$master_file) %>%
    dplyr::arrange(sample_id) %>% 
    dplyr::mutate(sample_id = as.character(sample_id)) %>%
    dplyr::mutate(
        diagnostic_group_coded = 
        ifelse(diagnostic_group == "CONTROL", 0, 
        ifelse(diagnostic_group == "MCI", 0.5, 1)))


#Load data with a bit of munging. Column names are better standardized on the full dataset
blood.df <- read.csv(args$blood_file) %>% 
    dplyr::arrange(sample_id) %>% 
    dplyr::mutate(sample_id = as.character(sample_id))

PCs.df <- read.csv(file.path(args$idir, paste0(args$chr, ".PCs.csv"))) %>% 
    dplyr::arrange(sample_id) %>% 
    dplyr::mutate(sample_id = as.character(sample_id))

#Pull phenotypes, PCs, and blood into one design matrix
#Subset to only have sample_ids also contained in M/Cov
df <- dplyr::inner_join(blood.df, master.df, by = "sample_id") %>%
        dplyr::inner_join(PCs.df, by = "sample_id")

#M(ethylated) reads and Cov(erage) from sequencing
M <- read_and_arrange("M")
Cov <- read_and_arrange("Cov")

# Correct sample IDs in the correct order
ordered_sample_ids <- c("chromStart", intersect(df$sample_id, names(M)))

# Don't get rid of the distinct
df <- df %>% 
    dplyr::filter(sample_id %in% ordered_sample_ids) %>% 
    dplyr::distinct()

# Rearrange columns to match data frame
M <- M %>% dplyr::select(all_of(ordered_sample_ids))
Cov <- Cov %>% dplyr::select(all_of(ordered_sample_ids))

#Check that ordering of samples in M/Cov is the same as in samplesheet
#should be in ascending order
if (all(df$sample_id == names(M)[-1])){
    print("Sample order is correct")
} else {
    warning("Sample order wrong!")
}

median_coverage <- apply(X=Cov[ ,-"chromStart"], FUN=median, MARGIN=1)
# min_coverage <- apply(X=Cov[ ,-"chromStart"], FUN=min, MARGIN=1)

keepix <- (median_coverage >= args$coverage_threshold)
pp <- sum(keepix) / length(keepix)
print(paste0(round(100*pp, 2), "% of sites meet coverage theshold"))

# Filter M and Cov based on median coverage threshold
M.mat <- as.matrix(M[keepix, -"chromStart"])
Cov.mat <- as.matrix(Cov[keepix, -"chromStart"])

# create bs seq object, needs chromosome identiifer, methylated reads, and unmethylated reads
bs <- BSseq(chr = rep(args$chr, nrow(M.mat)), 
            pos = M$chromStart[keepix],
            M = M.mat, 
            Cov =Cov.mat, 
            sampleNames = names(M)[-1])

# Not neccessary for running serial, but when parallel
# it helps to clean up
invisible(gc())

save(bs, df, file = ofile)
print(paste0("Wrote out ", ofile))

#END