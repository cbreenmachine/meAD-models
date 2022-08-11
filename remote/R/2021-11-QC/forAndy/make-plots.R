library(tidyverse)
library(data.table)
library(wiscR)

ifile <- "../../../data/batch01/pool01-group02/04-extract/110.tsv"
ifile.2 <- "./110.cpg.cov"

# Read and compute coverage informative for methylation
DT <- fread(ifile)
DT$coverage <- DT$methylated + DT$unmethylated

# Filtering not done in BCF
DT <- DT[2 < coverage & coverage < 80, methylation:=methylated / coverage]

p <- DT %>%
    ggplot(aes(x = methylation)) +
    geom_histogram(color = "black", bins = 50, aes(y = ..density..)) +
    geom_density(color = "red", adjust = 3, size = 1.6) +
    wiscR::light_theme() +
    xlim(c(0,1)) +
    xlab("Methylation estimate") +
    ylab("Density") +
    ggtitle("Density of methylation estimates")

wiscR::save_plot(p, "methylation-density.png")



#--> Coverage
DT <- read.table(ifile.2)

p <- DT %>%
    ggplot(aes(x = V1)) +
    geom_histogram(color = "black", bins = 25, aes(y = ..density..)) +
    geom_density(color = "red", adjust = 10, size = 1.6) +
    wiscR::light_theme() +
    xlab("Coverage") +
    ylab("Density") +
    xlim(c(0, 50)) +
    ggtitle("Density of coverage at CpG sites")

wiscR::save_plot(p, "coverage-density.png")
