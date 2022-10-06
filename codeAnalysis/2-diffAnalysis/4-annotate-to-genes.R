
suppressPackageStartupMessages({
    library(viridis)
    library(tidyverse)
    library(data.table)
    library(argparse)
    library(cowplot)
    library(ggsci)
    library(ggbio)
    library(EnsDb.Hsapiens.v86)
}) 

parser <- ArgumentParser(description='')
parser$add_argument('--idir', default = "../../data/07-counts-by-chrom-imputed-subset2/", help = "input directory.")
parser$add_argument('--LFDR', default = 0.05, help = "LFDR cutoff")
parser$add_argument('--upstream', default = 2000, help = "")
args <- parser$parse_args()


#TODO: fix out of range warning
expand_genes <- function(gr, upstream = 2000, downstream = 300) {
# https://bioinformatics.stackexchange.com/questions/4390/expand-granges-object-different-amounts-upstream-vs-downstream
  strand_is_minus = strand(gr) == "-"
  on_plus = which(!strand_is_minus)
  on_minus = which(strand_is_minus)
  start(gr)[on_plus] = start(gr)[on_plus] - upstream
  start(gr)[on_minus] = start(gr)[on_minus] - downstream
  end(gr)[on_plus] = end(gr)[on_plus] + downstream
  end(gr)[on_minus] = end(gr)[on_minus] + upstream
  return(gr)
}

#--> Load data and create Granges object
df <- fread(file.path(args$idir, "DMPs.bed"), 
            select = c("chrom", "chromStart", "lfdr", "LOAD.minus.control"))
gr <- df %>% dplyr::filter(lfdr < args$LFDR) %>%
        dplyr::mutate(chrom = str_remove(chrom, "chr")) %>%
        dplyr::rename(start=chromStart) %>% 
        dplyr::mutate(end=start+1) %>% 
        GRanges()


#--> Ensembl database and filtered for coding genes only
ensdb <- EnsDb.Hsapiens.v86
ensdb.subset <- genes(ensdb, filter = GeneBiotypeFilter('protein_coding'))
ensdb.expanded <- expand_genes(ensdb.subset)


#--> Intersect
overlaps <-findOverlaps(ensdb.expanded, gr)
hits <- queryHits(overlaps)

# Min and max should be in [0, 20,000]
summary(hits)

top_genes <- mcols(ensdb.expanded[hits]) %>%
    as.data.frame() %>%
    group_by(gene_name) %>%
    summarize(Count = n()) %>%
    arrange(-Count) %>%
    head(10)







#--> Top of plot, plots the gene introns/exons
plot_gene_model <- function(ensdb, gene_of_interest){
  
  p <- autoplot(ensdb, ~ symbol == gene_of_interest, 
                stat="reduce", label=F, 
                color="brown", fill="brown") +
    theme_void() +
    ggtitle(paste0(gene_of_interest, " gene model")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab("") +
    scale_x_continuous(position = "top") 
  
  return(p@ggplot)
}


#--> Bottom of plot, scatter of 
plot_points <- function(df, chrom, left, right){
  p <- df %>%
    dplyr::filter(chrom == chrom) %>%
    dplyr::filter(chromStart > left & chromStart + 1 < right) %>%
    dplyr::mutate(Significance = ifelse(lfdr < 0.05, "LFDR < 0.05", "LFDR > 0.05")) %>%
    ggplot(aes(x = chromStart, y = -log10(lfdr) * sign(LOAD.minus.control), color = Significance)) +
    geom_point(size = 0.5) +
    theme_minimal() +
    xlab("Genomic position") +
    scale_color_nejm() +
    xlim(c(left, right)) +
    ylab(expression('-log'[10]*'(LFDR)')) +
    theme(strip.text.x = element_blank(),
          strip.background = element_rect(colour="white", fill="white"),
          legend.position=c(.9,.87)
    ) 
  
  return(p)
}


#--> Driver function
make_figure <- function(ensdb, df, gene){
  print(gene)
  subset <- select(ensdb, GeneNameFilter(gene)) %>%
            dplyr::filter(SEQNAME == .$SEQNAME[1])
  
  left <- min(subset$TXSEQSTART)
  right <- max(subset$TXSEQEND)
  chrom <- paste0("chr", subset$SEQNAME[1])
  
  p1 <- plot_gene_model(ensdb, gene)
  p2 <- plot_points(df, chrom, left, right)
  p <-  cowplot::plot_grid(p1, p2, ncol = 1, axis = "b", align = "v") +
        theme(plot.background = element_rect(fill = "white", colour = "white"))
  
  return(p)
}

for (ii in 1:nrow(top_genes)) {
    gene <- top_genes[ii, 'gene_name']
    count <- top_genes[ii, 'Count']

    tryCatch({
    p <- make_figure(ensdb, df, gene)

    ofile <- file.path(args$idir, "figs", "genes-with-DMPs", paste0(gene, "-", count, ".png"))
    dir.create(dirname(ofile), showWarnings = F)
    cowplot::save_plot(filename = ofile, p)
    })
}



gene <- "BRCA1"
count <- 0
p <- make_figure(ensdb, df, gene)

ofile <- file.path(args$idir, "figs", "genes-with-DMPs", paste0(gene, "-", count, ".png"))
dir.create(dirname(ofile), showWarnings = F)
cowplot::save_plot(filename = ofile, p)
