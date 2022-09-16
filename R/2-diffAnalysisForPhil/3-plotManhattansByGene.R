# plotManhattans

library(BiocManager)
options(repos = BiocManager::repositories())
library(data.table)
library(tidyverse)
library(ggbio)
library(EnsDb.Hsapiens.v86)
library(cowplot)
library(ggsci)
library(scales)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--ifile", default= "../../dataDerived/methylBedImputedSubset1/withoutMCI/tests/has_risk_allele/DMPs.bed", help='Input DMPs file')
args <- parser$parse_args()

#Read data, need DMPs and ensembl database with coordinates
DMPs.df <- fread(args$ifile)  %>% dplyr::mutate(chrom = factor(chrom))
ensdb <- EnsDb.Hsapiens.v86 # DO NOT CHANGE
genes.df <- read_csv("../../dataReference/app.genes.csv")


#Style stuff--font size, colors
global_size = 22

colors <- pal_jama("default", alpha=1)(4)
c_gene <- colors[4]
c_sig <- colors[2]
c_not_sig <- colors[3]


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



#--> Driver function
make_figure <- function(igene, df=genes.df, save=FALSE){
  row <- df %>% dplyr::filter(gene == igene)

  chr <- dplyr::pull(row, 'chrom')
  start <- dplyr::pull(row, 'start')
  end <- dplyr::pull(row, 'end')
  strand <- dplyr::pull(row, 'strand')
  
  #TODO: add left/right sliders
  p1 <- plot_gene_model(chr, igene, strand=strand)

  autoplot_range <- layer_scales(p1)$x$range$range

  p2 <- plot_points(chr, start, end) + xlim(autoplot_range)
  p <- join_via_cow(p1, p2)
  
  if (save){
    cowplot::save_plot(filename=paste0(igene, ".png"), plot = p, 
                       base_width = 13, base_height = 6)
  }
  
  return(p)
}

genes <- c("PSEN1", "PSEN2", "APP", "APOE")

for (g in genes){
  make_figure(g, save=TRUE)
}

