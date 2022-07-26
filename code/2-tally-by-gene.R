library(tidyverse)
library(data.table)
library(argparse)
library(cowplot)
library(ggsci)
library(ggbio)
library(EnsDb.Hsapiens.v86)


parser <- ArgumentParser()
parser$add_argument('--ifiles', 
                    default = "data/DMPs-1.bed,data/DMPs-2.bed")
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
# Dummy set with just the cpg loci. Set lfdr to 0 to pass filter
cpgs.df <- fread(all_files[1], select = c("chrom", "chromStart")) %>% dplyr::mutate(lfdr = 0)

df.1 <- fread(all_files[1], select = c("chrom", "chromStart", "lfdr", "LOAD.minus.control"))
df.2 <- fread(all_files[2], select = c("chrom", "chromStart", "lfdr", "LOAD.minus.control"))


#--> First tally number of CpGs per gene
cpgs.gr <- df_to_gr(cpgs.df)
cpgs.overlap <-findOverlaps(ensdb.expanded, cpgs.gr)
cpgs.hits <- queryHits(cpgs.overlap)

#--> Now for the subsets
gr.1 <- df_to_gr(df.1)
rm(df.1)
gr.2 <- df_to_gr(df.2)
rm(df.2)
gc()


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

cpgs.per.gene <- pull_genes(cpgs.hits) %>% dplyr::rename(n_cpg = Count)
genes.1 <- pull_genes(hits.1) %>% dplyr::mutate(subset="subset1")
genes.2 <- pull_genes(hits.2) %>% dplyr::mutate(subset="subset2")

gene_count.df <- full_join(genes.1, genes.2, by="gene_name") %>%
  full_join(cpgs.per.gene, by = "gene_name") %>%
  dplyr::transmute(gene_name, 
                   n_sig_subset1 = Count.x,
                   n_sig_subset2 = Count.y,
                   n_cpg) %>%
  dplyr::mutate(percent_sig_subset1 = n_sig_subset1 / n_cpg,
                percent_sig_subset2 = n_sig_subset2 / n_cpg) %>%
  dplyr::mutate(percent_sig_avg = 0.5 * (percent_sig_subset1 + percent_sig_subset2))

gene_count.df %>%
  dplyr::filter(n_cpg > 100) %>%
  arrange(-percent_sig_avg) %>%
  head(100)


p <- gene_count.df %>%
  ggplot(aes(x = n_sig_subset1, y = n_sig_subset2, label=gene_name)) +
  geom_text(size = 2) +
  theme_bw()

cowplot::save_plot(p, filename="gene_counts.png")

