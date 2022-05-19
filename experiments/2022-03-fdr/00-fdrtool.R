library(tidyverse)
library(fdrtool)

ifile <- "../../analysis/2022-01-diff-analysis/chr18/smooth-150-PCs-2/models.RData"
load(ifile)

df <- data.frame(chr = test.cohort$chr, pos = test.cohort$pos, 
                stat = test.cohort$stat, pvals = test.cohort$pvals) %>% drop_na()


png("pval-before-fdr.png")
hist(df$pvals)
dev.off()

#--> Before
png("hist-of-stats.png")
hist(df$stat)
dev.off()

# After normalization
zz <- (df$stat - mean(df$stat)) / sd(df$stat)
png("hist-of-normalized.png")
hist(zz)
dev.off()

#--> FDTR Tool
out <- fdrtool(zz, statistic = "normal", plot = FALSE, cutoff.method = "locfdr")

png("qval.png")
hist(out$qval)
dev.off()


png("fdr.png")
hist(out$lfdr)
dev.off()


png("pval-after-fdr.png")
hist(out$pval)
dev.off()
