library(data.table)
library(tidyverse)
library(ggbio)
library(EnsDb.Hsapiens.v86)
library(cowplot)
library(ggsci)
library(scales)
library(argparse)
library(magrittr)

#CONSTANTS
sub_cols <- c("chrom", "chromStart", "lfdr", "LOAD.minus.control")
alpha <- 0.05
log10.alpha <- log10(alpha)

# COLORS
global_size = 24
colors <- pal_jama("default", alpha=1)(7)
c_gene <- colors[4] # Color of gene model

c_hyper <- colors[2] 
c_not_sig <- colors[7]
c_hypo <- colors[3] 



read_and_process <- function(ifile){
  fread(ifile, select=sub_cols) %>% 
    dplyr::mutate(chrom = factor(chrom),
                  position = chromStart,
                  signed.logged.lfdr = sign(LOAD.minus.control) * log10(lfdr),
                  status = ifelse(signed.logged.lfdr < log10.alpha, "Hypo", 
                               ifelse(signed.logged.lfdr > -log10.alpha, "Hyper", "NotSig"))) %>%
    return()
  
}

DMP.1 <- read_and_process("../data/DMPs-1.bed")
DMP.2 <- read_and_process("../data/DMPs-2.bed")
genes.df <- fread("../data/app.genes.csv") %>% dplyr::mutate(gene = factor(gene), chrom = factor(chrom))
ensdb <- EnsDb.Hsapiens.v86 # DO NOT CHANGE

DMPs.df <- full_join(DMP.1, DMP.2, by=c("chrom", "position"))

rm(DMP.1, DMP.2)
gc()


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


global_size = 22
#--> Bottom of plot, scatter of 
sig_levels = c("Hyper-Hyper", "Hyper-NotSig", "Hyper-Hypo",
               "NotSig-Hyper", "NotSig-NotSig", "NotSig-Hypo",
               "Hypo-Hyper", "Hypo-NotSig", "Hypo-Hypo")
# sig_colors = c("")

plot_points <- function(chrom, start, end, df=DMPs.df){
  p <- df %>%
    drop_na() %>%
    dplyr::filter(chrom == chrom) %>%
    dplyr::filter(position > start & position + 1 < end) %>%
    unite('Significance', c('status.x', 'status.y'), sep = "-", remove = FALSE) %>%
    dplyr::mutate(Significance = factor(Significance, levels=sig_levels),
                  SigSubset2 = factor(status.y, levels=c("Hyper", "NotSig", "Hypo"))) %>%
    ggplot(aes(x = position, y = signed.logged.lfdr.x, color = SigSubset2)) +
    geom_point(size = 3, alpha=0.8) +
    xlab("Genomic position") +
    xlim(c(start, end)) +
    ylim(c(-3, 3)) +
    geom_hline(yintercept = 0, color = "black", alpha = 0.5) +
    geom_hline(yintercept = log10(0.05), color = "black", alpha = 0.3) +
    geom_hline(yintercept = -log10(0.05), color = "black", alpha = 0.3) +
    ylab(expression('-log'[10]*'(LFDR)')) +
    scale_color_manual(values = c(c_hyper, c_not_sig, c_hypo)) +
    theme_dmp() %>%
    return(p)
}


join_via_cow <- function(p1, p2){
  p <- cowplot::plot_grid(p1, p2, ncol = 1, axis = "b",
                          align = "v", rel_heights = c(1, 2)) +
    theme(plot.background = element_rect(fill = "white", colour = "white"))
  return(p)
}



#--> Driver function
make_figure <- function(ichrom=NULL, igene, left=NULL,right=NULL, df=genes.df, save=FALSE){
  
  if (!is.null(ichrom)){
    row <- df %>% dplyr::filter(chrom==ichrom, gene == igene)
  } else {
    row <- df %>% dplyr::filter(gene == igene)
    ichrom <- dplyr::pull(row, 'chrom')
  }
  

  start <- dplyr::pull(row, 'start')
  end <- dplyr::pull(row, 'end')
  strand <- dplyr::pull(row, 'strand')
  
  #TODO: add left/right sliders
  p1 <- plot_gene_model(ichrom, igene, strand=strand)
  p2 <- plot_points(ichrom, start, end)
  
  if (!is.null(left) & !is.null(right)){
    print("Custom x lims")
    p1 <- p1 + xlim(c(left, right))
    p2 <- p2 + xlim(c(left, right))
  }
  p <- join_via_cow(p1, p2)
  
  if (save){
    cowplot::save_plot(filename=paste0("../figs/comparison-chicagos/comp-", igene, ".png"), plot = p, 
                       base_width = 9, base_height = 6)
  }
  return(p)
}


gc()
top_100 <- fread("../data/DMPs-by-gene-comparison-top100.csv")
all_genes <- c("RYR1", "RYR2", "RYR3", "APP", "PSEN1", "PSEN2", "APOE", top_100$gene_name)

for (gg in all_genes){
  
  try(make_figure(igene=gg, save = TRUE))
  Sys.sleep(time = 1)
}



#END

