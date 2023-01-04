# 0_prepare_DSS_chrom.R
# takes M and Cov matrices, alongside phenotypes and outputs an R data store with bs object (input for DSS)
# and design matrix (df)

suppressPackageStartupMessages({
    library(data.table)
    library(argparse)
    library(DSS)
    library(magrittr)
    library(bsseq)
})

parser <- ArgumentParser()
parser$add_argument("--idir", default= "../DataRaw/MCovRawSplit/", help='Directory to run DSS on')
parser$add_argument("--odir", default= "../DataDerived/ControlLOAD/Inputs/", help='Directory to run DSS on')
parser$add_argument("--blood_file", default="../DataRaw/methylBedForDXM/bloodCellPropsDXM.csv")
parser$add_argument("--master_file", default="../DataRaw/masterSamplesheet.csv")
parser$add_argument("--chr", default= "chr22", help='Chromosome to run DSS on')
# parser$add_argument("--keep_mci", help="TODO: NOT IMPLEMENTED")
parser$add_argument("--median_threshold", default=5, help="Keep sites with at least () median coverage")
parser$add_argument("--min_threshold", default=0, help="Keep sites with at least () minimum coverage")

args <- parser$parse_args()

#TODO--send this code to generate master... don't want to have to repeat this
master.df <- read.csv(args$master_file) %>%
    dplyr::arrange(sample_id) %>% 
    dplyr::mutate(sample_id = as.character(sample_id)) %>%
    dplyr::filter(diagnostic_group != "MCI") %>%
    dplyr::mutate(
        diagnostic_group_coded = 
        ifelse(diagnostic_group == "CONTROL", 0, 1))

# Check that the samples in the directory have phenotype and vice versa
idir.samples <- 
    unique(stringr::str_split_fixed(list.files(args$idir), "\\.", 3)[ ,2])

# All samples that are in the directory and sampleshseet
valid.samples <- intersect(master.df$sample_id, idir.samples)
load.samples <- intersect(
    master.df$sample_id[master.df$diagnostic_group == "LOAD"], 
    idir.samples)

ctrl.samples <- intersect(
    master.df$sample_id[master.df$diagnostic_group == "CONTROL"], 
    idir.samples
    )


read_wrapper <- function(s){
    #Creates input name like ../Data/chr22.100.bed
    # and reads it
    dt <- fread(file.path(args$idir, paste0(args$chr, ".", s ,".bed")))
    dt$sample <- as.numeric(s) # numeric quiets warning
    dt
}

data <- do.call(rbind, lapply(X=valid.samples, FUN=read_wrapper))

pivot_me <- function(data, value) {
    # Gets into Sample \times Position matrix of 
    # M or Cov
    # "value" is "methylated" or "coverage"
    keepcols <- c("chromStart", "sample", value)

    data %>% 
        dplyr::select(all_of(keepcols)) %>%
        tidyr::pivot_wider(values_from = value, names_from = sample) %>%
        tibble::column_to_rownames("chromStart")
}

M <- pivot_me(data, "methylated")
Cov <- pivot_me(data, "coverage")

# Set to zero so we can compute summary stats
Cov.zeroed <- Cov
Cov.zeroed[is.na(Cov.zeroed)] <- 0

min.coverage <- apply(Cov.zeroed, FUN=min, MARGIN=1)
med.coverage <- apply(Cov.zeroed, FUN=median, MARGIN=1)

# Pass both filters...
keepix <- ((min.coverage >= args$min_threshold) & 
           (med.coverage >= args$median_threshold))

# No need to impute if there's a minimum filter
M.filt <- M[keepix, ]
Cov.filt <- Cov[keepix, ]

# Number of samples
# N.load <- length(load.samples)
# N.ctrl <- length(ctrl.samples)
N <- ncol(M.filt)

create_filler_mask <- function(DT, ctrl.samples, load.samples, N){
    # N is number of columns (samples)
    row.means.load <- floor(rowMeans(DT[ , load.samples], na.rm = T))
    row.means.ctrl <- floor(rowMeans(DT[ , ctrl.samples], na.rm = T))

    # Expand to make it a matrix we can mask
    filler <- DT
    
    filler[ , load.samples] <- row.means.load
    filler[ , ctrl.samples] <- row.means.ctrl

    filler
}

# Check that this worked
# unique(unlist(filler[1000, ]))

# LOAD/Control row means expanded to be the same shape as M 
M.filler <- create_filler_mask(M.filt, ctrl.samples, load.samples, N)
Cov.filler <- create_filler_mask(Cov.filt, ctrl.samples, load.samples, N)

# Masks for indexing
mask <- is.na(Cov.filt)

# Fill withe means
M.filt[mask] <- M.filler[mask]
Cov.filt[mask] <- Cov.filler[mask]

# sum(is.na(Cov.filt))

invisible(gc())

# Derived the output directory from the name of the input
odir <- args$odir
dir.create(odir, showWarn = F, recursive = T)
ofile <- file.path(odir, paste0("input-", args$chr, ".RData"))


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
# after processing
M <- M.filt
Cov <- Cov.filt

# Correct sample IDs in the correct order
ordered_sample_ids <- c(intersect(df$sample_id, names(M)))

# Don't get rid of the distinct
df <- df %>% 
    dplyr::filter(sample_id %in% ordered_sample_ids) %>% 
    dplyr::distinct()

# Rearrange columns to match data frame
M <- M %>% dplyr::select(all_of(ordered_sample_ids))
Cov <- Cov %>% dplyr::select(all_of(ordered_sample_ids))

#Check that ordering of samples in M/Cov is the same as in samplesheet
#should be in ascending order
if (all(df$sample_id == names(M))){
    print("Sample order is correct")
} else {
    warning("Sample order wrong!")
}


# create bs seq object, needs chromosome identiifer, methylated reads, and unmethylated reads
bs <- BSseq(chr = rep(args$chr, nrow(M)), 
            pos = as.numeric(rownames(M)),
            M = as.matrix(M), 
            Cov = as.matrix(Cov), 
            sampleNames = names(M))

# Order positions out of abundance of caution
bs <- orderBSseq(bs)

# Not neccessary for running serial, but when parallel
# it helps to clean up
invisible(gc())

save(bs, df, file = ofile)
print(paste0("Wrote out ", ofile))

#END