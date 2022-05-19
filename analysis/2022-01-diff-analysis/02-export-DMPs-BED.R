# 02-call-DMRs.R
# Runs Dispersion Shrinkage Estimation method on methylation data (one chromosome at a time)
suppressPackageStartupMessages({
    library(data.table)
    library(tidyverse)
    library(argparse)
    library(fdrtool)
    library(VennDiagram)
    library(RColorBrewer)
})

# Argument parsing -----------------------------------------------------------
parser <- ArgumentParser()
parser$add_argument("--idir", default= "./chr18/smooth-150-PCs-2/", help="directory with models.RData")
parser$add_argument("--alpha", default= 0.01, help="local false discovery rate cutoff")
args <- parser$parse_args()


# Load data
load(file.path(args$idir, "models.RData"))

# Infer chromosome...
CHR <- test.cohort$chr[1]

# Coerce to normal data.frame
df <- test.cohort
class(df) <- "data.frame"
df <- drop_na(df) # a few hundred-couple thousand rows dropped

# given the processed test.cohort
zz <- (df$stat - mean(df$stat)) / sd(df$stat)
out <- fdrtool(zz, statistic = "normal", plot = FALSE, cutoff.method = "locfdr")

# Overwride pvals with corrected ones
df$pvals <- out$pval

# The rest are additional columns
df$pval <- out$pval
df$lfdr <- out$lfdr 
df$qval <- out$qval


###################### OUTPUT BED FILE ############################

ofile <- file.path(args$idir, "DMPs.bed")

df %>%
    filter(lfdr < args$alpha) %>%
    transmute(chr, start = pos - 1, end = pos + 1, lfdr = lfdr) %>%
    write_tsv(ofile)







################# PLOT VENN DIAGRAM ##############################

low_q <- df$pos[df$qval < alpha]
low_lfdr <- df$pos[df$lfdr < alpha]


venn.diagram(
        x = list(low_q, low_lfdr),
        category.names = c("Q" , "FDR"),
        filename = 'venn_diagramm.png',
        output=TRUE,
        
        # Output features
        imagetype="png" ,
        height = 480 , 
        width = 480 , 
        resolution = 300,
        compression = "lzw",

        col=c("#440154ff", '#21908dff'),
        fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
                
        # Numbers
        cex = .6,
        fontface = "bold",
        fontfamily = "sans",
        
)