suppressPackageStartupMessages({
    library(data.table)
    library(tidyverse)
    library(argparse)
    library(GenomicRanges)
    library(parallel)
    library(EnsDb.Hsapiens.v86)
    library(clusterProfiler)
    library(org.Hs.eg.db)
})

parser <- ArgumentParser()
parser$add_argument("--regions_file", default= "../../dataSummaries/controlLOAD-pc2-ct5/dmrs.DSS.bed", help='Directory to run DSS on')
parser$add_argument("--pvals_file", default= "../../dataSummaries/controlLOAD-pc2-ct5/pvals.bed", help='Directory to run DSS on')
parser$add_argument("--diff_file", default="../../dataSummaries/loadMinusControl.bed")
parser$add_argument("--odir", default="../../figs/controlLOAD-pc2-ct5/")
args <- parser$parse_args()


.expand_genes <- function(gr, upstream = 2000, downstream = 300) {
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

.pad_chrom <- function(s="chr1"){
    # chr1 --> chr01
    t <- str_remove(s, "chr") %>% 
        str_pad(width = 2, pad = "0") 
    paste0("chr", t) %>% return()
}

pvals.df <- fread(args$pvals_file)
regions.df <- fread(args$regions_file)
diff.df <- fread(args$diff_file) %>% dplyr::mutate(chrom = .pad_chrom(chrom))


# odir <- file.path(args$odir, "thresholdedDMRs")
# dir.create(odir, showWarn=F, recurs=T)


# Type conversion to get access to intersection...
.convert <- function(df){
    # data frame --> GRanges
    makeGRangesFromDataFrame(
        df,
        keep.extra.columns=T,
        starts.in.df.are.0based = T)
}

regions.gr <- .convert(regions.df)
diff.gr <- .convert(diff.df)
pvals.gr <- .convert(pvals.df)

# Exapmle on one
ix <- 1
sub.pvals.gr <- subsetByOverlaps(pvals.gr, regions.gr[ix])
sub.diff.gr <- subsetByOverlaps(diff.gr, regions.gr[ix])

tmp.df <- cbind(data.frame(sub.pvals.gr), 
                LOAD.minus.CONTROL = sub.diff.gr$LOAD.minus.CONTROL) %>% 
    dplyr::mutate(p.corrected = -log10(p.corrected)) %>%
    pivot_longer(cols = c("p.corrected", "LOAD.minus.CONTROL", "diagnostic_group_coded"), 
        names_to = "stat.type", 
        values_to = "stat")

p <- tmp.df %>%
    ggplot(aes(x = start, y = stat)) +
    geom_point() +
    geom_smooth(se=F) +
    facet_wrap(.~stat.type, scales="free", ncol=1) +
    xlab("Genomic position")

cowplot::save_plot("tmp.png", p) 



# Subset to DMPs in regions, so we don't have to query with 
# 25 million eveyr time
# diff.sub.gr <- subsetByOverlaps(diff.gr, regions.gr)


# Prints the LOAD.minus.CONTORL for CpGs in first entry of regions.gr
subsetByOverlaps(diff.sub.gr, regions.gr[1, ])$LOAD.minus.CONTROL

get_mean_in_range <- function(i){
    cpgs <- subsetByOverlaps(diff.sub.gr, regions.gr[i, ])
    return(mean(cpgs$LOAD.minus.CONTROL))
}

means <- mclapply(X=1:length(regions.gr), 
                FUN=get_mean_in_range, 
                mc.cores=12)

means.vec <- as.vector(do.call(rbind, means))
regions.gr$means <- means.vec



# Plot these DMRs
# for (ii in 1:length(regions.top.gr)){
#     my_range <- ranges(regions.top.gr[ii, ])
#     my_start <- start(my_range)
#     my_end <- end(my_range)
#     my_chrom <- as.character(seqnames(regions.top.gr[ii, ]))

#     p <- diff.df %>%
#         dplyr::filter(chrom==my_chrom) %>%
#         dplyr::filter(chromStart < my_end, chromStart > my_start) %>%
#         ggplot(aes(x = chromStart, y = LOAD.minus.CONTROL)) +
#         geom_smooth(se=F) +
#         geom_point() +
#         theme_light() +
#         xlab("Position") +
#         ggtitle(paste0("Mean diff: ", regions.top.gr[ii, ]$means)) +
#         ylim(c(-1,1))

#     ofile <- file.path(odir, paste0( "groupMean-", my_chrom, "-", my_start, ".png") )
#     cowplot::save_plot(file=ofile, p)
#     print(ii)
# }



# Subset regions to 
cut <- 0.025
regions.top.gr <- regions.gr[abs(regions.gr$means) > cut, ]

#--> Ensembl database and filtered for coding genes only
ensdb <- EnsDb.Hsapiens.v86
ensdb.subset <- genes(ensdb, filter = GeneBiotypeFilter('protein_coding'))
ensdb.subset <- ensdb.subset[seqnames(ensdb.subset) %in% 1:22]

seqlevelsStyle(ensdb.subset) <- "dbSNP"
seqlevelsStyle(regions.top.gr) <- "dbSNP"

ensdb.expanded <- .expand_genes(ensdb.subset, upstream=2000, downstream=300)

overlaps <-findOverlaps(ensdb.expanded, regions.top.gr)
hits <- queryHits(overlaps)

top_genes <- unique(as.vector(ensdb.expanded[hits, ]$gene_id))

top_gene_names <- unique(as.vector(ensdb.expanded[hits, ]$gene_name))
write.table(sort(top_gene_names), "2022-10-03-ADPD-topGenes.txt", quote=F, row=F, col=F)

go.out <- clusterProfiler::enrichGO(top_genes, keyType="ENSEMBL",
            OrgDb = 'org.Hs.eg.db', ont="ALL")


go.df <- as.data.frame(go.out) %>%
    arrange(-p.adjust)

write_csv(go.df, "2022-10-03-ADPD-geneOntology.csv")

sink("2022-10-03-summary.txt")
paste("Number of regions: ", length(regions.top.gr))
paste("Cutoff in mean differences: ", cut)

sink()