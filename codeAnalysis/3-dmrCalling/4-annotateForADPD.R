suppressPackageStartupMessages({
    library(data.table)
    library(tidyverse)
    library(argparse)
    library(ggbio)
    library(EnsDb.Hsapiens.v86)
    library(cowplot)
    library(ggsci)
    library(scales)
    
})

ensdb <- EnsDb.Hsapiens.v86 # DO NOT CHANGE
genes.df <- read_csv("../../dataReference/app.genes.csv")



parser <- ArgumentParser()
parser$add_argument("--regions_file", default= "../../dataSummaries/DMPs-controlLOAD/combp.regions-p.bed", help='Directory to run DSS on')
parser$add_argument("--dmps_file", default= "../../dataSummaries/DMPs-controlLOAD/DMPs.bed", help='Directory to run DSS on')
parser$add_argument("--diff_file", default= "../../dataSummaries/loadMinusControl.bed", help='Directory to run DSS on')
parser$add_argument("--odir", default= "../../outputsForAbstracts/2023-ADPD", help='output directory')
args <- parser$parse_args()

dir.create(args$odir, recurs=T, showWar=F)
genes.df <- read_csv("../../dataReference/app.genes.csv")


#Style stuff--font size, colors
global_size = 22

colors <- pal_jama("default", alpha=1)(4)
c_gene <- colors[4]
c_sig <- colors[2]
c_not_sig <- colors[3]



my_subset <- function(data, cc, left, right){
    dplyr::filter(data, chrom==cc) %>%
    dplyr::filter(chromStart >= left && chromEnd <= right) %>%
    return()
}

regions.df <- read_table(args$regions_file) %>% 
    dplyr::filter(n_probes > 15) %>%
    dplyr::rename(chrom = `#chrom`, chromStart = start, chromEnd = end)

diff.df <- fread(args$diff_file)
dmps.df <- fread(args$dmps_file)

my_subset(dmps.df, "chr6", 0, 1000000)



# Theme -------------------------------------------------------------------
theme_dmp = function(){
  theme_minimal() %+replace%
  theme(strip.text.x = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position=c(.9,.87),
        text = element_text(size=global_size),
        plot.title = element_text(hjust = 0.5)
  ) 
}


# Gene model function -----------------------------------------------------
plot_gene_model <- function(ichrom, igene, strand='*', db=ensdb){
  # Top of plot, plots the gene introns/exons
  iseq <- str_remove(ichrom, pattern='chr0|chr')
  my_filter <- AnnotationFilter(~ symbol == igene & seq_name == iseq)
  
  p <- autoplot(db, my_filter, 
                stat="reduce", label=F, 
                color= c_gene, fill= c_gene) +
    ggtitle(paste0(igene, " gene model (", strand, " strand)")) +
    ylab("") +
    scale_x_continuous(position = "top") +
    theme_dmp()
    
  return(p@ggplot)
}


# Not, is
sig_levels <- c("LFDR \u2265 0.05", "LFDR < 0.05")

#--> Bottom of plot, scatter of 
plot_points <- function(chr, start, end, df=DMPs.df){
  df.tmp <- df %>%
    dplyr::filter(chrom == chr) %>%
    dplyr::filter(chromStart > start & chromStart + 1 < end) %>%
    dplyr::mutate(signed.logged.p = - sign(has_risk_allele) * log10(lfdrs))
  
  p <- df.tmp %>%
    dplyr::mutate(Significance = 
                    factor(
                      ifelse(abs(signed.logged.p) <  -log10(0.05),
                           sig_levels[1], sig_levels[2]),
                      levels = sig_levels)) %>%
    ggplot(aes(x = chromStart, y = signed.logged.p, color = Significance)) +
    geom_point(size = 3) +
    xlab("Genomic position") +
    scale_color_manual(values = c(c_sig, c_not_sig)) +
    xlim(c(start - 500, end + 500)) +
    ylim(c(-3, 3)) +
    geom_hline(yintercept = 0, color = "black", alpha = 0.5) +
    geom_hline(yintercept = log10(0.05), color = "black", alpha = 0.3) +
    geom_hline(yintercept = -log10(0.05), color = "black", alpha = 0.3) +
    ylab(expression('-log'[10]*'(LFDR)')) +
    theme_dmp()

  return(p)
}


join_via_cow <- function(p1, p2){
  p <- cowplot::plot_grid(p1, p2, ncol = 1, axis = "b",
                          align = "v", rel_heights = c(1, 2)) +
    theme(plot.background = element_rect(fill = "white", colour = "white"))
  return(p)
}



