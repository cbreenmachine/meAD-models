# 3-exportDMPsAsBED.R
# 
suppressPackageStartupMessages({
    library(data.table)
    library(tidyverse)
    library(argparse)
    library(ComplexHeatmap)
    library(viridis)
    library(circlize)
    library(ggsci)
})

parser <- ArgumentParser()
parser$add_argument("--idir", default= "../../dataDerived/analysis-controlLOAD/test-diagnostic-group-coded/", help='Where are the models stored')
parser$add_argument("--odir", default= "../../figs/2022-paper/heatmaps/", help='Where are the figs saved')
args <- parser$parse_args()

dir.create(args$odir, showWarn=F)

# Gets us the covariates in the right order
#TODO: change this!!!! May not always be the same order when we run other covariates
load("../../dataDerived/analysis-controlLOAD/test-diagnostic-group-coded/output.chr22.RData")

#Load and filter DMRs
dmrs <- fread(file.path(args$idir, "dmrs.all.bed"))
dmrs.pass <- dmrs %>% 
    dplyr::filter(dist_to_nearest_gene == 0)

subset_pi_mat <- function(df, dmr.start, dmr.end){
    # This way of filtering ensures same number of CpGs
    df %>% 
        dplyr::filter(
            start >= dmr.start,
            start <= dmr.end
        ) %>% 
        dplyr::select(-c("chr", "start", "end")) %>%
        data.matrix() %>% 
        t()
}

make_plot_title <- function(my.chr, dmr.start, dmr.end, my.gene, n.cpgs, n.dmps){
    paste0("DMR in ", my.gene,
                " (", my.chr, ":", 
                format(dmr.start, big.mark=",", sci=F), "-",
                format(dmr.end, big.mark=",", sci=F), ")\n",
                n.dmps, " out of ", n.cpgs, " CpGs are DMPs")

}


make_annotation <- function(load.cat, sex.cat){
    HeatmapAnnotation(
        DG = load.cat, 
        Sex = sex.cat,
        which = "row", 
        col = list(
            DG = c("C"="#868686FF", "L"="#8843F2"),
            Sex = c("M" = "#7AA6DCFF", "F"="#CD534CFF"))
        )
}



plot_hm_matrix <- function(sub.pi, my.title, ha){

    n.cpgs <- ncol(sub.pi)

    if (n.cpgs > 10){
        n.cpgs <- 1
    }

    ComplexHeatmap::Heatmap(
        sub.pi, col = magma(50), 
        column_order = 1:ncol(sub.pi), #turn off clustering on CpGs
        heatmap_legend_param = list(title = "pi", at = c(0, 0.5, 1)),
        left_annotation = ha, column_title = my.title,
        row_km = 2, # KNN for rows (samples)
        column_split = 1:n.cpgs,
        border = F)
}

for (my.chr in unique(dmrs$chr)){

    chr.dmrs <- dplyr::filter(dmrs.pass, chr == my.chr)
    if (nrow(chr.dmrs) == 0){next}

    pi.df <- fread(file.path(args$idir,  "experimentSummary", paste0("pi.", my.chr, ".bed")))

    for (i in 1:nrow(chr.dmrs)){
        # Extract info for title
        dmr.start <- chr.dmrs$start[i]
        dmr.end <- chr.dmrs$end[i]
        my.gene <- chr.dmrs$nearest_gene_name[i]
    
        # Subset data
        sub.pi <- subset_pi_mat(pi.df, dmr.start, dmr.end)
         
        # tally
        n.cpgs <- ncol(sub.pi)
        n.dmps <- chr.dmrs$n.sig[i]
        
        my.title <- make_plot_title(my.chr, dmr.start, dmr.end, my.gene, n.cpgs, n.dmps)
        
        # Test to make sure the order is correct
        all(rownames(design.df) == rownames(sub.pi))
        load.cat <- ifelse(design.df$diagnostic_group == 0, "C", "L")
        sex.cat <- ifelse(design.df$sex == "Male", "M", "F")

        # Kill rownames to de clutter heatmap
        rownames(sub.pi) <- NULL
        ha <- make_annotation(load.cat, sex.cat)

        # Save plot
        ofile <- file.path(args$odir, paste0(my.gene, "-", my.chr, ":", dmr.start, "-", dmr.end, "-heatmap.png"))
        print(ofile)

        png(ofile, width = 8, height = 8, unit="in", res = 150)
        print(plot_hm_matrix(sub.pi, my.title, ha))
        dev.off()
    }
}

#END