library(tidyverse)
library(wiscR)
library(ggsci)
library(argparse)
library(viridis)

study_file_path <- "../../data/study-data/experimental-design.tsv"
study.df <- read.table(study_file_path, sep="", header=TRUE)

parser <- ArgumentParser(description='')
parser$add_argument('--idir', default = "../../data/07-counts-by-chrom-imputed-subset1/", help = "input directory.")
parser$add_argument('--chrom', default = "chrom1", help = "input directory.")
args <- parser$parse_args()

meta.file <- list.files(args$idir, pattern="master", full=TRUE)[1]

# Make output and collect chromosomal files
samples.df <- read_csv(meta.file)

# Data Munging, need the PCs and colors to be in same df
df <- file.path(args$idir, "chr1.PCs.bed") %>% read_tsv()
df.2 <- merge(df, samples.df, by.x="sample_id", by.y="RA_id") %>% 
    mutate(sample = sample_id) %>%
    left_join(study.df, by = 'sample') 

cols.to.factor <- c("apoe_e1", "apoe_e2", "education", "excluded_from_R01", "dxgrp_bin", "age_at_visit")
df.2[ cols.to.factor] <- lapply(df.2[cols.to.factor], factor)


plot_pca <- function(df, color_by_str, odir, chrom, label = FALSE){
    # 2d scatterplot of first two principal components
    # colored by 'color_by_str' which should be string
   print(color_by_str)
    p <- df %>%
        drop_na() %>%
        ggplot(aes_string(x = "PC1", y = "PC2", color = color_by_str)) +
        geom_point(size = 6, alpha = 0.8) +
        wiscR::light_theme() +
        ggtitle("First two PCs") 

    num_uniq <- length(unique(df[[color_by_str]]))
    if (is.numeric(df[[color_by_str]]) & num_uniq > 2) {
        p <- p + scale_fill_viridis()
        
    } else if (length(unique(df[[color_by_str]])) < 10) {
        p <- p + scale_color_npg()
    
    }

    if (label){
        # adds '101', '432', etc. to the points
        # font sizes and colors are a bit wacky, but helps identify outliers
        p <- p + geom_text(aes(label=sample),hjust=0, vjust=0, size = 16)
        # prepend <- paste0(prepend, "_labeled")
    }
    
    ofile <- file.path(odir, "figs", paste0(chrom, ".PCs.", color_by_str, ".png"))
    wiscR::save_plot(p, ofile)
}


plot_wrapper <- function(color_by_str){
    plot_pca(df.2, color_by_str, args$idir, args$chrom)
}

variables <- c("diagnostic_group", "race_primary", "apoe_e1", "apoe_e2", 
    "education", "excluded_from_R01", "Charlson.Co.Morbidity.Index.Score", "dxgrp_bin",
    "age_at_visit", "hispanic", "gender", "machine", "pool", "group", "batch")
lapply(variables, plot_wrapper)

# plot_pc_by_name <- function(ff){
#     # grab the beginning of the important part of the string for naming
#     #chr <- str_replace(basename(ff), ".csv", "")
#     df <- read_tsv(ff, col_types = cols() ) %>% left_join(samples.df, by = "sample")

#     all_vars <- names(df)
#     good_vars <- all_vars[!str_detect(all_vars, "PC|sample")]

#     for (vv in good_vars){
#         plot_pca(df, vv)
#         print(vv)
#     }
    
#     plot_pca(df, "mean_methylation", label=TRUE)
#     print(ff)
# }

# lapply(all_files, plot_pc_by_name)


# # plot_scree <- function(ff){

# #     ve.df <- read_csv(ff)
# #     chr <- str_remove(basename(ff), ".csv")

# #     p <- ve.df %>%
# #         dplyr::filter(chrom == chr) %>%
# #         arrange(PC) %>%
# #         mutate(cum_var_explained = cumsum(var_explained)) %>%
# #         ggplot(aes(x = PC, y = cum_var_explained)) +
# #         geom_point(size = 7) +
# #         geom_line(size = 1.3) +
# #         ylim(c(0, 0.5)) +
# #         xlab("Number of PCs") +
# #         ylab("Cumulative var explained (%)") +
# #         ggtitle(chr) +
# #         wiscR::light_theme()

# #     wiscR::save_plot(p, file.path(args$odir, paste0("scree-", chr, ".png")))
# # }


