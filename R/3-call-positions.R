# 3-exportDMPsAsBED.R
suppressPackageStartupMessages({
    library(data.table)
    library(magrittr)
    library(argparse)
})

parser <- ArgumentParser()
parser$add_argument("--ifile", default= "../DataDerived/ControlLOAD/test-diagnostic-group-coded/Summaries/pvals.bed", help='Where are the models stored')
parser$add_argument("--ofile", default= "../DataDerived/ControlLOAD/test-diagnostic-group-coded/Summaries/DMPs.bed", help='Where are the models stored')
parser$add_argument("--ofile_dmps", default= "../DataDerived/ControlLOAD/test-diagnostic-group-coded/Summaries/DMPs-high-effect.bed", help='Where are the models stored')
parser$add_argument("--fdr_cut", default = 0.05)
parser$add_argument("--effect_cut", default = 0.1)
args <- parser$parse_args()

LFDR.CUT <- as.numeric(args$fdr_cut)
EFFECT.CUT <- as.numeric(args$effect_cut)


ofile.1 <- args$ofile

a <- stringr::str_pad(LFDR.CUT * 100, 2, pad="0")
b <- EFFECT.CUT * 100
# ofile.2 <- file.path(args$odir, paste0("dmps-lfdr", a,"-effect",b, ".bed"))
ofile.2 <- args$ofile_dmps  

df <- fread(args$ifile)

dmps.all <- df %>% dplyr::filter(lfdr < LFDR.CUT)
dmps.filt <- dmps.all %>% dplyr::filter(abs(pi.diff) >= EFFECT.CUT)

fwrite(dmps.all, ofile.1, sep="\t")
fwrite(dmps.filt, ofile.2, sep="\t")
# END