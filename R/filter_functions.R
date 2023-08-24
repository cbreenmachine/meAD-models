

.read_bed <- function(idir, chr, sample){
    #Creates input name like Data/chr22.100.bed
    # and reads it
    dt <- fread(file.path(idir, paste0(chr, ".", sample ,".bed")))
    dt$sample <- as.numeric(sample) # numeric quiets warning
    dt
}

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


read_and_pivot_routine <- function(idir, chr, samples){
    # Read and stack all sample-specific data files
    .read_bed_wrapper <- function(sample){
        .read_bed(idir, chr, sample)
    }
    data <- do.call(rbind, lapply(X = samples, FUN = .read_bed_wrapper))

    M <- pivot_me(data, "methylated")
    Cov <- pivot_me(data, "coverage")

    return(list(M = M, Cov = Cov))
}

# # Set to zero so we can compute summary stats

filter_and_impute <- function(m_cov){
    m <- m_cov[["M"]]
    cov <- m_cov[["Cov"]]

    cov_zeroed <- data.table::copy(cov) 
    cov_zeroed[is.na(cov_zeroed)] <- 0

    # We'll use the number of samples as a divisor
    n_samples <- ncol(cov_zeroed)
    
    percent_zero <- rowSums(Cov.zeroed == 0) / n_samples
    percent_nonzero <- 1 - percent_zero

    # median coverage calculation happens independently 
    # from zero
    med_cov <- apply(cov_zeroed, FUN=median, MARGIN=1)

# # Pass both filters...
# # Note that a CpG can pass if (e.g.) more than half of samples
# # have positive coverage, and the median is above 5. But there could
# # still be samples with zero coverage.
# # Consider 10, 10, 10, 10, 10, 0, 0.
# keepix <- ((percent.nonzero >= args$percent_nonzero_threshold) & 
#            (med.coverage >= args$median_threshold))

# percent.keep <- round(100 * sum(keepix) / length(keepix), 2)
# paste0("Keeping ", percent.keep, "%")

# # Save the M and Cov matrices if requested
# # These have the missing values; useful to check
# # manually for one or two chromosomes
# if (args$save_pivots){
#     ofile <- file.path(odir, paste0("pivot-", args$chr, ".RData"))
#     print(paste0("Saving ", ofile))

#     save(M, Cov, percent.keep, file = ofile)
# }

# # No need to impute if there's a minimum filter
# M.filt <- data.table::copy(M)[keepix, ]
# Cov.filt <- data.table::copy(Cov)[keepix, ]

}

# # Number of samples
# N <- ncol(M.filt)

# create_filler_mask <- function(DT, ctrl.samples, load.samples, N){
#     # Given the groups of samples, compute the mean for columns in 
#     # "load.samples" and then make a data table that looks like DT,
#     # but with columns imputed with mean

#     # N is number of columns (samples)
#     row.means.load <- round(rowMeans(DT[ , load.samples], na.rm = T))
#     row.means.ctrl <- round(rowMeans(DT[ , ctrl.samples], na.rm = T))

#     # Expand to make it a matrix we can mask
#     filler <- DT
    
#     filler[ , load.samples] <- row.means.load
#     filler[ , ctrl.samples] <- row.means.ctrl

#     filler
# }

# # LOAD/Control row means expanded to be the same shape as M 
# M.filler <- create_filler_mask(M.filt, ctrl.samples, load.samples, N)
# Cov.filler <- create_filler_mask(Cov.filt, ctrl.samples, load.samples, N)

# # sum(Cov.filler)
# # [1] 0 # as expected

# # Masks for indexing
# # mask <- is.na(Cov.filt)
# mask <- (Cov.filt == 0 | is.na(Cov.filt))
# percent.of.kept.imputed <- round(100 * sum(mask) / (nrow(M.filt) * ncol(M.filt)), 2)

# # Fill withe means
# M.filt[mask] <- M.filler[mask]
# Cov.filt[mask] <- Cov.filler[mask]

# # sum(Cov.filt == 0)
# # sum(is.na(Cov.filt))
# invisible(gc())

# #########################################
# ###### Compute Principal Components #####
# #########################################

# P <- t(M.filt / Cov.filt)
# dim(P)

# P.vars <- apply(P, 2, stats::var)

# # Get a numeric cut from the percentile
# cut <- quantile(P.vars, 1 - args$top_percent_most_variable)

# # Subet methylation matrix  
# P.top <- P[ ,(P.vars >= cut)]
# dim(P.top)

# # Compute PCs
# # P.top should be (and is) wide
# pca.out <- gmodels::fast.prcomp(P.top, scale. = T, center = T)
# PC.df <- as.data.frame(pca.out$x) 
# PC.df$sample_id <- rownames(PC.df)


# pc.ofile <- file.path(pc.odir, paste0("pcs-", args$chr, ".png"))
# png(pc.ofile)
# plot(PC.df$PC1, PC.df$PC2,
#     xlab = "PC1", ylab = "PC2")
# dev.off()


# #########################################
# ########## Samplesheet munging ##########
# #########################################

# #Load data with a bit of munging. Column names are better standardized on the full dataset
# blood.df <- read.csv(args$blood_file) %>% 
#     dplyr::arrange(sample_id) %>% 
#     dplyr::mutate(sample_id = as.character(sample_id))

# #Pull phenotypes, PCs, and blood into one design matrix
# #Subset to only have sample_ids also contained in M/Cov
# df <- dplyr::inner_join(blood.df, master.df, by = "sample_id") %>%
#         dplyr::inner_join(PC.df, by = "sample_id")

# #M(ethylated) reads and Cov(erage) from sequencing
# # after processing
# M <- M.filt
# Cov <- Cov.filt

# # Correct sample IDs in the correct order
# ordered_sample_ids <- c(intersect(df$sample_id, names(M)))

# # Don't get rid of the distinct
# df <- df %>% 
#     dplyr::filter(sample_id %in% ordered_sample_ids) %>% 
#     dplyr::distinct()

# # Rearrange columns to match data frame
# M <- M %>% dplyr::select(all_of(ordered_sample_ids))
# Cov <- Cov %>% dplyr::select(all_of(ordered_sample_ids))

# #Check that ordering of samples in M/Cov is the same as in samplesheet
# #should be in ascending order
# if (all(df$sample_id == names(M))){
#     print("Sample order is correct")
# } else {
#     warning("Sample order wrong!")
# }

# # create bs seq object, needs chromosome identiifer, methylated reads, and unmethylated reads
# bs <- BSseq(chr = rep(args$chr, nrow(M)), 
#             pos = as.numeric(rownames(M)),
#             M = as.matrix(M), 
#             Cov = as.matrix(Cov), 
#             sampleNames = names(M))

# # Order positions out of abundance of caution
# bs <- orderBSseq(bs)

# # Not neccessary for running serial, but when parallel
# # it helps to clean up
# invisible(gc())

# # Derived the output directory from the name of the input
# save(bs, df, pca.out, percent.keep, percent.of.kept.imputed, med.coverage, file = ofile)
# print(paste0("Wrote out ", ofile))
# #END


extract_filtered_data <- function(bs){
    
}