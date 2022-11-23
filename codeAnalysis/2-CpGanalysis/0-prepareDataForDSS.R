# 0-prepareDataForDSS.R
# takes M and Cov matrices, alongside phenotypes and outputs an R data store with bs object (input for DSS)
# and design matrix (df)

suppressPackageStartupMessages({
    library(data.table)
    library(tidyverse)
    library(argparse)
    library(DSS)
})

parser <- ArgumentParser()
parser$add_argument("--idir", default= "../../dataDerived/methylBedImputed-controlMCILOAD/", help='Directory to run DSS on')
parser$add_argument("--blood_file", default="../../dataDerived/methylBedForDXM/bloodCellPropsDXM.csv")
parser$add_argument("--master_file", default="../../dataDerived/masterSamplesheet.csv")
parser$add_argument("--chr", default= "chr22", help='Chromosome to run DSS on')

args <- parser$parse_args()

# Derived the output directory from the name of the input
odir <- str_replace(args$idir, "methylBedImputed", "analysis")
dir.create(odir, showWarn = F, recursive = T)
ofile <- file.path(odir, paste0("input.", args$chr, ".RData"))

if (!file.exists(ofile)){

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

    #Pull phenotypes, PCs, and blood into one design matrix
    #Subset to only have sample_ids also contained in M/Cov
    df <- inner_join(blood.df, master.df, by = "sample_id") %>%
            inner_join(PCs, by = "sample_id")

    #M(ethylated) reads and Cov(erage) from sequencing
    M <- fread(file.path(args$idir, paste0(args$chr, ".M.bed")))
    Cov <- fread(file.path(args$idir, paste0(args$chr, ".Cov.bed")))

    # Correct sample IDs in the correct order
    ordered_sample_ids <- c("chromStart", intersect(df$sample_id, names(M)))

    # Don't get rid of the distinct
    df <- df %>% dplyr::filter(sample_id %in% ordered_sample_ids) %>% distinct()

    # Rearrange columns to match data frame
    M <- M %>% dplyr::select(all_of(ordered_sample_ids))
    Cov <- Cov %>% dplyr::select(all_of(ordered_sample_ids))

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

    gc()

    if (! exists("bs")){stop(paste0(ofile, " could not be created"))}

    save(bs, df, file = ofile)
    print(paste0("Wrote out ", ofile))

} else {print(paste0(ofile, " already exists, delete if you want to re-run!"))}