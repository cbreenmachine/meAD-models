library(tidyverse)
library(wiscR)
library(ggsci)
library(argparse)
library(viridis)

# study_file_path <- "../../data/study-data/experimental-design.tsv"
# study.df <- read.table(study_file_path, sep="", header=TRUE)

parser <- ArgumentParser(description='')
parser$add_argument('--idir', default = "../../data/07-counts-by-chrom-imputed-subset2/", help = "input directory.")
parser$add_argument('--chr', default = "chr1", help = "which chromosome to run on")
args <- parser$parse_args()

# meta.file <- list.files(args$idir, pattern="master", full=TRUE)[1]
master.df <- read_csv(file.path(args$idir, "master.csv"))
blood.df <- read_csv(file.path(args$idir, "blood-cell-composition.csv"))
pc.df <- read_tsv(file.path(args$idir, paste0(args$chr, ".PCs.bed")))

# Data Munging, need the PCs and colors to be in same df
df <- merge(master.df, blood.df, by.x="RA_id", by.y="sample") %>% 
    merge(pc.df, by.x = 'RA_id', by.y = "sample_id") %>%
    select(-c(study_id, dxgrp_bin))

cols.to.factor <- c("apoe_e1", "apoe_e2", "education", "excluded_from_R01", "sex", "hispanic")
df[cols.to.factor] <- lapply(df[cols.to.factor], factor)


plot_pca <- function(df, color_by_str, odir, chr, label = FALSE){
    # 2d scatterplot of first two principal components
    # colored by 'color_by_str' which should be string
   print(color_by_str)
    p <- df %>%
        drop_na() %>%
        ggplot(aes_string(x = "PC1", y = "PC2", color = color_by_str)) +
        geom_point(size = 6, alpha = 0.8) +
        wiscR::light_theme() +
        ggtitle("First two PCs") 

    if (is.numeric(df[[color_by_str]])) {
        p <- p + scale_color_viridis(option = "magma")
        
    } else {
        p <- p + scale_color_npg()
    
    }

    if (label){
        # adds '101', '432', etc. to the points
        # font sizes and colors are a bit wacky, but helps identify outliers
        p <- p + geom_text(aes(label=sample),hjust=0, vjust=0, size = 16)
        # prepend <- paste0(prepend, "_labeled")
    }
    
    dir.create(file.path(odir, "principal-components"))
    ofile <- file.path(odir, "figs/principal-components",
            paste0(chr, ".PCs.", color_by_str, ".png"))
    wiscR::save_plot(p, ofile)
}


plot_wrapper <- function(color_by_str){
    plot_pca(df, color_by_str, args$idir, args$chr)
}

# variables <- c("diagnostic_group", "race_primary", "apoe_e1", "apoe_e2", 
#     "education", "excluded_from_R01", "Charlson.Co.Morbidity.Index.Score", "dxgrp_bin",
#     "age_at_visit", "hispanic", "gender", "machine", "pool", "group", "batch")
variables <- names(df)[-1]
lapply(variables, plot_wrapper)


p <- df %>%
        drop_na() %>%
        ggplot(aes_string(x = "PC1", y = "PC2", label="RA_id")) +
        geom_text(size = 10) +
        wiscR::light_theme() +
        ggtitle("First two PCs")

ofile <- file.path(args$idir, "figs", paste0(args$chr, ".PCs.text.png"))
wiscR::save_plot(p, ofile)