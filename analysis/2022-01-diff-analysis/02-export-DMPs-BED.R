# 02-call-DMRs.R
# Runs Dispersion Shrinkage Estimation method on methylation data (one chromosome at a time)
suppressPackageStartupMessages({
    library(data.table)
    library(tidyverse)
    library(fdrtool)
})


all.files <- system("find . -name models.RData | grep smooth-150-PCs-2", intern = T)

.correct_pvals <- function(df){
    zz <- (df$stat - mean(df$stat)) / sd(df$stat)
    out <- fdrtool(zz, statistic = "normal", plot = FALSE, cutoff.method = "locfdr")

    df$pvals.adj <- out$pval
    df$lfdr <- out$lfdr 
    df$qval <- out$qval
    return(df)
}


.load_data <- function(path){
    # Load from path like "./chr8/smooth-150-PCs-2/models.RData"
    load(path)
    df <- data.frame(chr = test.cohort$chr, pos = test.cohort$pos, 
                stat = test.cohort$stat) %>% drop_na() 
    return(df)
}

# .wrapper <- function(path){ .load_data(path) %>% .correct_pvals() %>% return()}
.wrapper <- function(path){ .load_data(path) %>% return()}
df <- do.call("rbind", lapply(all.files, .wrapper))  %>% .correct_pvals() %>% mutate(p = pvals.adj)

df.2 <- df %>% 
    dplyr::mutate(chrom = chr, chromStart = pos - 1, chromEnd = chromStart + 1, direction = sign(stat)) %>% 
    dplyr::select(-pos) %>%
    dplyr::select(chrom, chromStart, chromEnd, lfdr, p, qval, direction)


filter_and_write <- function(df, s, cut=0.05){
    
    sig <- paste0("0", as.character(round(cut * 100)))
    ofile <- paste0("DMPs.", s, sig, ".bed")
    print(ofile)

    df %>% filter(get(s) < cut) %>%
        fwrite(ofile, sep = "\t")
}

for (cut in c(0.01, 0.05)){

    for (s in c("lfdr", "p", "qval")){
        filter_and_write(df.2, s, cut)
    }

}

ofile <- "DMPs.bed"
fwrite(df.2, ofile, sep="\t")