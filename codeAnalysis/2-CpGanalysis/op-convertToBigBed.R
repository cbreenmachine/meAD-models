library(data.table)
library(rtracklayer)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

ifile <- "../../dataSummaries/controlLOAD-pc2-ct5/pvals.bed"
ofile <- "../../dataSummaries/tracksUCSC/2022-10-14-pvals.bb"

df <- fread(ifile)

# Assign groups, take mean over 10 CpGs?

df$score <- -log10(df$p.corrected) * sign(df$diagnostic_group_coded)
ix <- !is.na(df$score)

gr <- makeGRangesFromDataFrame(
    df[ix, c("chrom", "chromStart", "chromEnd", "score")],
    keep.extra = T, ignore.strand = T, starts.in.df.are.0based = T
    )

# seqlevels(gr) <- c(paste0("chr", 1:22), "chrX", "chrY")
seqinfo(gr) <- seqinfo(txdb)[seqnames(seqinfo(gr))]

export(gr, ofile)
