# 02-call-DMRs.R
# Runs Dispersion Shrinkage Estimation method on methylation data (one chromosome at a time)
suppressPackageStartupMessages({
    library(data.table)
    library(tidyverse)
    library(bacon)
    library(argparse)
    library(cowplot)
    library(DSS)
    library(fastqq)
    library(latex2exp)
})

parser <- ArgumentParser()
parser$add_argument("--idir", default= "../../dataDerived/models-controlMCILOAD-pc0-ct0/diagnostic_group_coded/", help='Where are the models stored')
parser$add_argument("--odir", default= "../../dataSummaries/controlMCILOAD-pc0-ct0/", help='directory to store DMPs')
parser$add_argument("--fig_dir", default= "../../figs/pvals-controlMCILOAD-pc0-ct0/", help='figure outputs')
parser$add_argument("--title", default= "2 PCs, 5 cell types", help='figure title')
args <- parser$parse_args()

# Output directories
dir.create(args$odir, recursive=T, showWarn=F)
dir.create(args$fig_dir, recursive=T, showWarn=F)
all.files <-list.files(args$idir, pattern = "*RData", full=T)

#Derive the degrees of freedom for p-value correction
load(all.files[1])
N <-nrow(df)
n_cov <- length(names(beta.df))
degrees_freedom <- N - n_cov + 1

# test.cohort has chr,pos,stat,pvals,fdrs as columns and is DSS-sepcific
.load_data <- function(path){
    # Load from path like "./chr8/smooth-150-PCs-2/models.RData"
    load(path)

    # experimental design df
    N_this_chrom <- nrow(df)

    # overwrite with covariates and p-values
    df <- cbind(test.result, beta.df) %>% drop_na() 

    

    if (N_this_chrom != N){
        warning(paste0("Chromosome ", path, " has wrong number of samples: ", N_this_chrom))
    }
    # This allows us to do testing with DSS framework
    # class(df)[2] <- "DMLtest.multiFactor"
    return(df)
}


# How to do this efficiently... Loop thru chromosomes??
.load_wrapper <- function(path){ .load_data(path) %>% return()}
df <- do.call("rbind", lapply(all.files, .load_wrapper)) 


theme_set(theme_minimal(base_size=12) +
    theme(plot.background = element_rect(color="white", fill="white"))
)

# https://github.com/haowulab/DSS/blob/55cb1a05738a2ed02717af50b7b52828bc6b508d/R/DML.multiFactor.R#L192
# ofile <- file.path(args$fig_dir, "2022-09-30-testStatsSquaredRaw-forSunduz.png")

summarize_with_plot <- function(chi2_stat, pvals){

    l <- median(chi2_stat) / qchisq(0.5, 1)
    p.thinned <- fastqq::drop_dense(x=runif(n=length(pvals)), y=pvals)

    caption <- str_replace(str_replace(dss.formula, "type1", "\ntype1"),"PC1", "\nPC1")[2]

    qq.plot <- p.thinned %>%
        dplyr::mutate(x = -log10(x), y = -log10(y)) %>%
        ggplot(aes(x, y)) +
        geom_point() + 
        geom_line(aes(x=x,y=x), color="red") +
        xlab("Theoretical") +
        ylab("Observed") +
        ggtitle(TeX(sprintf(r'($\lambda = %.2f$)', l))) +
        labs(caption = paste0("N samples =", N))
    
    hist.plot <- df %>%
        ggplot(aes(x = pvals)) +
        geom_histogram(bins=20, aes(y=..density..)) +
        xlab("Unadjusted p-value") +
        ylab("Density") +
        ggtitle(args$title) +
        labs(caption = caption)

    joined <- cowplot::plot_grid( hist.plot, qq.plot)
}

chi2_stat <- df$stat ^ 2
pvals <- df$pvals

joined <- summarize_with_plot(chi2_stat, pvals)

ofile <- file.path(args$fig_dir, "2022-10-02-qqWithHist.png")
cowplot::save_plot(ofile, joined)

# z_scores <- (df$stat - mean(df$stat)) / sd(df$stat)
z_scores <- df$stat

#--> Bacon routine
bc <- bacon(z_scores, niter=20000L)

text_ofile <- file.path(args$odir, "bacon.txt")
sink(text_ofile)
print(bc)
sink()

chi2_stat <- qchisq(1-pval(bc),1)
chi2_stat <- chi2_stat[is.finite(chi2_stat)]

pvals <- pval(bc)

joined.bacon <- summarize_with_plot(chi2_stat, pvals)

ofile <- ofile <- file.path(args$fig_dir, "2022-10-02-qqWithHist-postBacon.png")
cowplot::save_plot(ofile, joined.bacon)


# ofile <- file.path(args$odir, "DMPs.bed")
# fwrite(DMPs.df, ofile, sep="\t")


DMPs.df <- df %>% 
    dplyr::transmute(chrom = chr, chromStart = pos - 1, chromEnd = pos, 
            diagnostic_group_coded, p.raw = pvals)

DMPs.df$p.adj <- pval(bc)
  
ofile <- file.path(args$odir, "DMPs.bed")
fwrite(DMPs.df, ofile, sep="\t")

# Outputs for comb-p
.pad_chrom <- function(s="chr1"){
    # chr1 --> chr01
    t <- str_remove(s, "chr") %>% 
        str_pad(width = 2, pad = "0") 
    paste0("chr", t) %>% return()
}



# combp wants slightly different format
# use adjusted p-values here
combp.df <- DMPs.df %>%
    transmute(chrom=.pad_chrom(chrom), start=chromStart, end=chromEnd, p=p.adj) %>%
    arrange(chrom, start)

ofile <- file.path(args$odir, "DMPs.forCombP.bed")
fwrite(combp.df, ofile, sep="\t")

#END