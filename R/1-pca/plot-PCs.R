library(tidyverse)
library(wiscR)
library(ggsci)
library(argparse)
library(viridis)

master.path <- "../../dataDerived/masterSamplesheet.csv"
master.df <- read_csv(master.path, show_col_types = FALSE)


cat.vars <- c("diagnostic_group", "machine", "pool", "sex", "race_primary", "hispanic",
                "batch", "source", "apoe_e1", "apoe_e2", "excluded_from_r01", "smoke")
num.vars <- c("bmi", "age_at_visit", "education","charlson_co_morbidity_index_score")

parser <- ArgumentParser(description='')
parser$add_argument('--idir', default = "../../dataDerived/methylBedImputed/withMCI/", help = "input directory.")
parser$add_argument('--odir', default = "../../figs/prinCompsWithMCI/", help = "input directory.")
parser$add_argument('--chrom', default = "chr18", help = "which chromosome to run on")
args <- parser$parse_args()

dir.create(args$odir, showWarn=F, recursive=T)

# Data Munging, need the PCs and colors to be in same df
data <- file.path(args$idir, paste0(args$chrom, ".PCs.csv")) %>% 
    read_csv(show_col_types=F) %>%
    left_join(master.df, by="sample_id") %>%
    mutate_at(cat.vars, as.character) %>%
    mutate_at(num.vars, as.numeric)


plot_numeric <- function(color_by_str, df=data, odir=args$odir, chrom=args$chrom){
    p <- df %>%
        ggplot(aes_string(x = "PC1", y = "PC2", color = color_by_str)) +
        geom_point(size = 6, alpha = 0.8) +
        wiscR::light_theme() +
        scale_color_viridis() +
        ggtitle("First two PCs") +
        labs(caption = paste0("Principal components calculated on ", chrom)) +
        theme(plot.caption = element_text(hjust = 0),
              plot.background = element_rect(fill="white"))
    
    ofile <- file.path(odir, paste0(chrom, ".PCs.", color_by_str, ".png"))
    wiscR::save_plot(p, ofile)
}

for (v in num.vars){plot_numeric(v)}

plot_categorical <- function(color_by_str, df=data, odir=args$odir, chrom=args$chrom){
    
    p <- df %>%
            mutate(color_by_str = color_by_str) %>%
            ggplot(aes_string(x = "PC1", y = "PC2", color = color_by_str)) +
            geom_point(size = 6, alpha = 0.8) +
            wiscR::light_theme() +
            ggtitle("First two PCs") +
            labs(caption = paste0("Principal components calculated on ", chrom)) +
            theme(plot.caption = element_text(hjust = 0),
                plot.background = element_rect(fill="white"))
    
    if (length(unique(df[[color_by_str]])) <=8 ) {
        p <- p + scale_color_nejm()
    } else {
        # Nothing
    }
    
    ofile <- file.path(odir, paste0(chrom, ".PCs.", color_by_str, ".png"))
    wiscR::save_plot(p, ofile)
    
}

for (v in cat.vars){print(v); plot_categorical(v)}

p <- data %>%
    ggplot(aes(x = PC1, y = PC2, color = diagnostic_group, label=sample_id)) +
    geom_text(size=9) +
    wiscR::light_theme() +
    scale_color_igv() +
    ggtitle("First two PCs") +
    labs(caption = paste0("Principal components calculated on ", args$chrom)) +
    theme(plot.caption = element_text(hjust = 0),
            plot.background = element_rect(fill="white"))

ofile <- file.path(args$odir, paste0(args$chrom, ".PCs.text.png"))
wiscR::save_plot(p, ofile)


# BOX Plots

