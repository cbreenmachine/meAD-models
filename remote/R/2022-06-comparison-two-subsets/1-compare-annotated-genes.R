

library(viridis)
library(tidyverse)
library(data.table)
library(argparse)
library(cowplot)
library(ggsci)
library(ggbio)
library(EnsDb.Hsapiens.v86)


parser <- ArgumentParser()
parser$add_argument('--ifiles', 
    default = "../../data/07-counts-by-chrom-imputed-subset2/DMPs.bed,../../data/07-counts-by-chrom-imputed-subset2/DMPs.bed")
parser$add_argument('--LFDR', default = 0.05, help = "LFDR cutoff")
parser$add_argument('--upstream', default = 2000, help = "")
parser$add_argument('--downstream', default = 300, help = "")
args <- parser$parse_args()

all_files <- unlist(str_split(args$ifiles, ","))

#TODO: fix out of range warning
expand_genes <- function(gr, upstream = args$upstream, downstream = 300) {
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



df_to_gr <- function(df){
    df %>% dplyr::filter(lfdr < args$LFDR) %>%
            dplyr::mutate(chrom = str_remove(chrom, "chr")) %>%
            dplyr::rename(start=chromStart) %>% 
            dplyr::mutate(end=start+1) %>% 
            GRanges() %>%
            return()

}

#--> Ensembl database and filtered for coding genes only
ensdb <- EnsDb.Hsapiens.v86
ensdb.subset <- genes(ensdb, filter = GeneBiotypeFilter('protein_coding'))
ensdb.expanded <- expand_genes(ensdb.subset)


#--> Load data and create Granges object
df.1 <- fread(all_files[1], select = c("chrom", "chromStart", "lfdr", "LOAD.minus.control"))
df.2 <- fread(all_files[2], select = c("chrom", "chromStart", "lfdr", "LOAD.minus.control"))

gr.1 <- df_to_gr(df.1)
gr.2 <- df_to_gr(df.2)




#--> Intersect
overlaps.1 <-findOverlaps(ensdb.expanded, gr.1)
hits.1 <- queryHits(overlaps.1)

#--> Intersect
overlaps.2 <-findOverlaps(ensdb.expanded, gr.2)
hits.2 <- queryHits(overlaps.2)


pull_genes <- function(hits){

    mcols(ensdb.expanded[hits]) %>%
        as.data.frame() %>%
        group_by(gene_name) %>%
        summarize(Count = n()) %>%
        return()
}

genes.1 <- pull_genes(hits.1) %>% dplyr::mutate(subset="subset1")
genes.2 <- pull_genes(hits.2) %>% dplyr::mutate(subset="subset2")

gene_count.df <- full_join(genes.1, genes.2, by="gene_name")

p <- gene_count.df %>%
    ggplot(aes(x = Count.x, y = Count.y, label=gene_name)) +
    geom_text(size = 2) +
    theme_bw()
    
cowplot::save_plot(p, filename="gene_counts.png")

