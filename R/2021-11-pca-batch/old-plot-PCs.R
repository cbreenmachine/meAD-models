library(tidyverse)
library(wiscR)
library(argparse)
library(viridis)


parser <- ArgumentParser(description='Process some integers')
parser$add_argument('--samples_files', default = "../../data/meta/phenos-cleaned.csv", help = "input phenotypes to color by")
parser$add_argument('--idir', default = "../../data/prin-comps/", help = "input directory.")
parser$add_argument('--odir', default = "./figs/", help='where to store output figures.')
args <- parser$parse_args()

# Make output and collect chromosomal files
samples.df <- read_csv(args$samples_file) %>% dplyr::mutate(pool_group = paste(pool, group, sep="-"))
dir.create(args$odir, showWarnings = FALSE)
my_files <- list.files(path = args$idir, pattern =  "chr*.csv", recursive = FALSE, full = TRUE)
my_files <- my_files[str_detect(my_files, "chr")]

plot_pca <- function(df, color_by_str, label = FALSE){
    # 2d scatterplot of first two principal components
    # colored by 'color_by_str' which should be string
    #TODO: clean label string (title case and _ to " ")
    #TODO: 3 breaks in legend--min, mean, max? Why are they inconsistent
    p <- df %>%
        drop_na() %>%
        ggplot(aes_string(x = "PC1", y = "PC2", color = color_by_str)) +
        geom_point(size = 6, alpha = 0.8) +
        wiscR::light_theme() +
        ggtitle("First two PCs of methylation matrix") 
        #labs(caption = paste0("N observations: ", n, "\t", "N missing: ", n_na, "\n", "N imputed w/ local mean: ", n_imputed))

    if ((length(unique(df[[color_by_str]])) > 10) & (is.numeric(df[[color_by_str]]))) {
        p <- p + scale_color_viridis(option = "magma")
    }

    if (label){
        # adds '101', '432', etc. to the points
        # font sizes and colors are a bit wacky, but helps identify outliers
        p <- p + geom_text(aes(label=sample),hjust=0, vjust=0, size = 16)
        prepend <- paste0(prepend, "_labeled")
    }

    ofile <- file.path(args$odir, paste0(prepend, "_", color_by_str, ".png"))
    wiscR::save_plot(p, ofile)
}


plot_pc_by_name <- function(ff){
    # grab the beginning of the important part of the string for naming
    #chr <- str_replace(basename(ff), ".csv", "")
    df <- read_tsv(ff, col_types = cols() ) %>% left_join(samples.df, by = "sample")

    all_vars <- names(df)
    good_vars <- all_vars[!str_detect(all_vars, "PC|sample")]

    for (vv in good_vars){
        plot_pca(df, vv)
        print(vv)
    }
    
    plot_pca(df, "mean_methylation", label=TRUE)
    print(ff)
}

lapply(all_files, plot_pc_by_name)




plot_scree <- function(ff){

    ve.df <- read_csv(ff)
    chr <- str_remove(basename(ff), ".csv")

    p <- ve.df %>%
        dplyr::filter(chrom == chr) %>%
        arrange(PC) %>%
        mutate(cum_var_explained = cumsum(var_explained)) %>%
        ggplot(aes(x = PC, y = cum_var_explained)) +
        geom_point(size = 7) +
        geom_line(size = 1.3) +
        ylim(c(0, 0.5)) +
        xlab("Number of PCs") +
        ylab("Cumulative var explained (%)") +
        ggtitle(chr) +
        wiscR::light_theme()

    wiscR::save_plot(p, file.path(args$odir, paste0("scree-", chr, ".png")))
}


