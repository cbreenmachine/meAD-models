# 4-compareSetsOfDMPs.R
# Two outputs : bed file of common DMPs among 5 variables
# control/LOAD, RAVLT, abeta_42/40, p tau abeta42, p tau
# Outputs an upSet figure

suppressPackageStartupMessages({
    library(viridis)
    library(tidyverse)
    library(data.table)
    library(cowplot)
    library(ggvenn)
    library(ComplexHeatmap)
}) 

root_dir <- "../../dataDerived/analysis-controlLOAD/test-"
lfdr.cut <- 0.05

test.vars <- c("diagnostic-group-coded", "ravlt-long",
      "a-beta-42-40-bin", "p-tau-abeta42-bin", "p-tau-bin")

all.files <- paste0(root_dir, test.vars, "/experimentSummary/pvals.bed")


load_and_filter <- function(ifile){
    df <- fread(ifile) %>% 
        dplyr::filter(lfdr < lfdr.cut) %>% 
        dplyr::mutate(locus = paste0(chr, ":", chromStart))

    # grab 'diagnostic_group-coded'
    cov <- str_split_fixed(str_split_fixed(ifile, "test-", 2)[[2]], "/", 2)[[1]]
    df$locus
}

all.dmps <- lapply(all.files, load_and_filter)
names(all.dmps) <- test.vars

m <- make_comb_mat(all.dmps, mode = "intersect")
m.sub <- m[comb_degree(m) %in% c(1, 2, 5)]

#RENAME vars


ht_opt(
    legend_title_gp = gpar(fontsize = 12, fontface = "bold"), 
    legend_labels_gp = gpar(fontsize = 12), 
    heatmap_column_names_gp = gpar(fontsize = 8),
    heatmap_column_title_gp = gpar(fontsize = 10),
    heatmap_row_title_gp = gpar(fontsize = 8)
)




p <- ComplexHeatmap::UpSet(
    m.sub, 
    comb_order = order(-comb_size(m.sub)),
    top_annotation = upset_top_annotation(
        m.sub,
        axis_param = list(at = c(0, 1e5, 2e5, 3e5),
            labels = c("0", "100k", "200k", "300k")),
        height = unit(4, "cm")
    ),
    right_annotation = upset_right_annotation(
        m,
        axis_param = list(at = c(0, 1e5, 2e5, 3e5),
            labels = c("0k", "100k", "200k", "300k"),
            labels_rot = 0),
        width = unit(4, "cm")
    )
)

###################################
############ OUTPUTS ##############
###################################
ofile <- "../../figs/summaryPlots/2022-11-07-upSetControlLOAD.png"
png(file= ofile, width=6, height=4, unit="in", res=300, pointsize=6); draw(p); dev.off()

ofile <- "../../dataSummaries/commonDMPs.bed"
smoking.dmps <- Reduce(intersect, all.dmps) 

smoke.df <- data.frame(
    chr = str_split_fixed(smoking.dmps, ":", 2)[ ,1],
    start = as.numeric(str_split_fixed(smoking.dmps, ":", 2)[ ,2])
) %>% dplyr::mutate(end = start + 2)

head(smoke.df)
fwrite(smoke.df, ofile, sep="\t")
# END



# odir <- "../../figs/summaryPlots/"
# ofile <- file.path(odir, paste0(Sys.Date(),"-roseComparison.png"))

# p <- ggvenn(
#         all_dmps, 
#         fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
#         stroke_size = 0.5, set_name_size = 3) 
#     # theme(panel.background = element_rect(fill="white", color="white"),
#     #       plot.background = element_rect(fill="white", color="white"))

# cowplot::save_plot(ofile, p)


# # Control LOAD discrete vs continuous
# ofile <- file.path(odir, paste0(Sys.Date(),"-roseComparison-controlLOAD.png"))

# p <- ggvenn(
#         all_dmps[1:2], 
#         fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
#         stroke_size = 0.5, set_name_size = 4) 

# cowplot::save_plot(ofile, p)



# # Control MCI LOAD discrete vs continuous
# ofile <- file.path(odir, paste0(Sys.Date(),"-roseComparison-controlMCILOAD.png"))

# p <- ggvenn(
#         all_dmps[3:4], 
#         fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
#         stroke_size = 0.5, set_name_size = 4) 

# cowplot::save_plot(ofile, p)


# # Discrete 
# ofile <- file.path(odir, paste0(Sys.Date(),"-roseComparison-discrete.png"))

# p <- ggvenn(
#         all_dmps[c(1,3)], 
#         fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
#         stroke_size = 0.5, set_name_size = 4) 

# cowplot::save_plot(ofile, p)


# # Continuouse 
# ofile <- file.path(odir, paste0(Sys.Date(),"-roseComparison-continuous.png"))

# p <- ggvenn(
#         all_dmps[c(2,4)], 
#         fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
#         stroke_size = 0.5, set_name_size = 4) 

# cowplot::save_plot(ofile, p)
