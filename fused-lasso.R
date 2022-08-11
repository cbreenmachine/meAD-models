library(genlasso)
library(fdrtool)
library(tidyverse)

load("wgbs-load/analysis/2022-01-diff-analysis/chr18/smooth-150-PCs-2/models.RData")

df <- test.cohort 
class(df) <- "data.frame"

df <- df %>% drop_na()

pos <- df$pos
zz <- (df$stat - mean(df$stat)) / sd(df$stat)

out <- fdrtool(zz, statistic = "normal", plot = FALSE, cutoff.method = "locfdr")


y <- -log10(out$pval)


start <- 655121
N <- 1000

pp <- pos[(start-N):(start+N)]

mod <- fusedlasso1d(y = y[(start-N):(start+N)], pos = pp)

beta <- mod$bls
summary(beta)

png("path.png")
plot(pp, beta)
lines(pp, beta, type = "l")
dev.off()

