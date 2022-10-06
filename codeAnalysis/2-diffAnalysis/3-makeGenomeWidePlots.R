# Make three plots and stores them in input-directory/figs/
# 1. Manhattan
# 2. Volcano
# 3. Waffle plot (what percent of DMPs in intron, exon, etc.)

suppressPackageStartupMessages({
    library(viridis)
    library(tidyverse)
    library(data.table)
    library(wiscR)
    library(argparse)
    library(genomation)
    library(methylKit)
    library(waffle)
    library(ggsci)
    library(latex2exp)
}) 


parser <- ArgumentParser(description='')
parser$add_argument('--ifile', default = "../../dataSummaries/controlLOAD/DMPs.bed", help = "input directory.")

parser$add_argument('--odir', default = "../../figs/controlLOAD/", help = "output directory.")
args <- parser$parse_args()

odir <- file.path(args$odir, "genomeWide/")

my_pal <- pal_locuszoom()
p_cut_thinning <- 0.0001

alpha_bonferroni <- 0.01 / nrow(df)

# Load data --------------------------------------------------------------------
df <- fread(args$ifile)

# N <- sum(df$lfdrs < args$LFDR)
# print(paste("Num significant:", N))
dir.create(odir, showWarn=F, recursive=T)

# Pvalues before and after
ofile <- file.path(odir, "pvalsBeforeAdjustment.png")
png(ofile)
hist(df$p.raw, main="Raw p-values ")
dev.off()

ofile <- file.path(odir, "pvalsAfterAdjustment.png")
png(ofile)
hist(df$p.adj, main="Adjusted p-values (bacon)")
dev.off()



# Plot functions -------------------------------------------------------------------
# Thanks to https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/
.pad_chrom <- function(s="chr1"){
    # chr1 --> chr01
    t <- str_remove(s, "chr") %>% 
        str_pad(width = 2, pad = "0") 
    paste0("chr", t) %>% return()
}

.adjust_pos <- function(df){
    data.cum <- df %>% 
        group_by(chrom) %>% 
        summarise(max_bp = max(chromStart)) %>% 
        mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>% 
        dplyr::select(chrom, bp_add)

    df %>% 
        inner_join(data.cum, by = "chrom") %>% 
        mutate(bp_cum = chromStart + bp_add) %>% 
        return()
}

.generate_axis_set <- function(df){
    df %>% group_by(chrom) %>% 
        summarize(center = mean(bp_cum)) %>% return()
}

sub.df <- df %>% mutate(chrom = .pad_chrom(chrom)) %>% dplyr::arrange(chrom, chromStart) %>% .adjust_pos() 
axis.set <- .generate_axis_set(sub.df) %>% mutate(chrom = .pad_chrom(chrom)) %>% dplyr::arrange(chrom)

# Thin non-significant points
ix <- which(sample(df$p.adj) > p_cut_thinning)
ix.to.ignore <- ix[1:floor(length(ix) * 0.997)]


p <- sub.df[!ix.to.ignore, ]  %>%
  ggplot(aes(x = bp_cum, color = as.factor(chrom), y = -log10(p.adj))) +
  geom_point(alpha = 0.8) +
  scale_x_continuous(label = axis.set$chrom, breaks = axis.set$center) +
  theme_minimal() +
  ylab("-log10(P)") +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 0.5),
    plot.caption = element_text(hjust = 0, face= "italic"),
    panel.background = element_rect(colour = NA),
    plot.background = element_rect(colour = NA)) +
  scale_color_manual(values = rep(c("grey", my_pal(1)[1]), 12)) +
  xlab("Genomic position")

ofile <- file.path(odir, "manhattanAdjustedPFromBacon.png")
cowplot::save_plot(plot = p, filename=ofile, 
                   base_width = 7, base_height = 5)


# Volcano Plot ------------------------------------------------------------
ix <- which(sample(df$p.adj) > 0.2)
ix.to.ignore <- ix[1:floor(length(ix) * 0.997)]

# Base plot--handle points, colors, etc.
p <- df[!ix.to.ignore, ] %>%
  mutate(Significance = ifelse(p.adj < alpha_bonferroni, "Significant", "Not Significant")) %>%
  ggplot(aes(
    x = diagnostic_group_coded, y = -log10(p.adj), 
    color = Significance)) +
  geom_point()+
  scale_color_manual(values = c(my_pal(3)[3], "grey")) +
  xlab(TeX("Effect size ($\\hat{\\beta}_{hasLOAD}$)")) +
  ylab(TeX("-\\log10(P)")) +
  xlim(c(-2.75, 2.75)) +
  # ylim(c(0, 3.75)) +
  theme_minimal() +
  ggtitle("Volcano plot of DMPs") +
  theme(legend.position=c(.87, .85),
        text = element_text(size=18),
        plot.title = element_text(hjust = 0),
        plot.background = element_rect(color = "white", fill = "white"),
        legend.background = element_rect(color="grey", fill=alpha("white", 0.4)))


ofile <- file.path(odir,  "volcanoAdjustedPFromBcaon.png")
cowplot::save_plot(plot = p, filename=ofile, 
                   base_width = 7, base_height = 5)



for (ii in 1:20){
  
}


# Parameters for annotation
# x_begin <- 1.2
# x_end <- 2.5
# yy <- 3
# text_size = 4.5

# p.annotated <- p + annotate("segment", x = x_begin, y = yy, xend = x_end, yend = yy, 
#            color="black",
#            arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
#   annotate("text", x = x_begin, size=text_size, hjust=0, y = yy+0.3, label = "Hypermethylated") +
#   annotate("segment", x = -x_begin, y = yy, xend = -x_end, yend = yy, 
#            color="black",
#            arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
#   annotate("text", x = -x_end, y = yy+0.3, size=text_size, hjust=0, label = "Hypomethylated") 

# ofile <- file.path(odir, paste0(Sys.Date(), "-volcano-LFDR.png"))
# cowplot::save_plot(plot = p.annotated, filename=ofile, 
#                    base_width = 7, base_height = 5)


# Pie chart / annotation -----------------------------------------
gr <- df %>% dplyr::filter(lfdrs < args$LFDR) %>%
        dplyr::rename(start=chromStart) %>% 
        dplyr::mutate(end=start+1) %>% 
        GRanges()

gene.obj <- readTranscriptFeatures("../../dataReference/gencode.hg38.bed")
diff.ann <- annotateWithGeneParts(gr, gene.obj)


p <- data.frame(Percent = round(diff.ann@precedence)) %>% 
  rownames_to_column("Location") %>%
  ggplot(aes(fill = Location, values = Percent)) +
    geom_waffle(n_rows = 10, n_cols = 10, size = 0.25, color="white") +
    scale_fill_manual(name = NULL, values = pal_locuszoom()(4)) +
  coord_equal() +
  ggtitle("Annotation of DMPs") +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", colour = "white"),
        plot.title = element_text(hjust = 0.5, size=20),
        legend.text=element_text(size=16),
        legend.position = "bottom",
        plot.margin = margin(t=5, r=5,b=5,l=5)) +
  guides(fill=guide_legend(ncol=2))


ofile <- file.path(odir, paste0(Sys.Date(), "-waffle-DMPs-annotation.png"))
cowplot::save_plot(plot = p, filename=ofile, 
                   base_width = 7, base_height = 5)



# Summary stats
tmp <- df %>%
  group_by(chrom) %>%
  summarise(N.CpG = n(), N.Sig.CpG = sum(lfdrs < 0.05)) %>%
  mutate(chrom  = .pad_chrom(chrom)) %>%
  pivot_longer(cols = starts_with("N"), values_to = "Number", names_to = "Type")


p <- tmp %>%
  ggplot(aes(fill=chrom, y=Number, x=Type)) + 
    geom_bar(position="fill", stat="identity") +
    theme_minimal() +
    theme(
        text = element_text(size=18),
        plot.title = element_text(hjust = 0),
        plot.background = element_rect(color = "white", fill = "white"),
        legend.background = element_rect(color="grey", fill=alpha("white", 0.4)))

  
ofile <- file.path(odir, paste0(Sys.Date(), "-barChart-sigCpGsByChrom.png"))
cowplot::save_plot(plot = p, filename=ofile, 
                   base_width = 7, base_height = 5)

#END