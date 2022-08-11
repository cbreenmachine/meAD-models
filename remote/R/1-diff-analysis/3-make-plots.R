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
}) 

parser <- ArgumentParser(description='')
parser$add_argument('--idir', default = "../../data/07-counts-by-chrom-imputed-subset2/", help = "input directory.")
parser$add_argument('--LFDR', default = 0.05, help = "LFDR cutoff")
args <- parser$parse_args()

# Load data --------------------------------------------------------------------
df <- fread(file.path(args$idir, "DMPs.bed"), 
            select = c("chrom", "chromStart", "lfdr", "LOAD.minus.control"))
odir <- file.path(args$idir, "figs/")





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

ix <- which(sample(df$lfdr) > 0.2)
ix.to.ignore <- ix[1:floor(length(ix) * 0.997)]


p <- sub.df[!ix.to.ignore, ]  %>%
  ggplot(aes(x = bp_cum, color = as.factor(chrom), y = -log10(lfdr))) +
  geom_point(alpha = 0.8) +
  scale_x_continuous(label = axis.set$chrom, breaks = axis.set$center) +
  theme_minimal() +
  ylab("-log10(LFDR)") +
  geom_hline(yintercept = -log10(args$LFDR)) +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 0.5),
    plot.caption = element_text(hjust = 0, face= "italic"),
    panel.background = element_rect(colour = NA),
    plot.background = element_rect(colour = NA)) +
  scale_color_manual(values = rep(c("grey", "#C5050C"), 12)) +
  xlab("Genomic position")


ofile <- file.path(odir, paste0(Sys.Date(), "-manhattan-LFDR.png"))
cowplot::save_plot(plot = p, filename=ofile, 
                   base_width = 7, base_height = 5)


# Volcano Plot ------------------------------------------------------------
ix <- which(sample(df$lfdr) > 0.2)
ix.to.ignore <- ix[1:floor(length(ix) * 0.997)]

x_begin <- 1.2
x_end <- 2.2
yy <- 2.5

p <- df[!ix.to.ignore, ] %>%
  mutate(Significance = ifelse(lfdr < args$LFDR, "LFDR < 0.05", "LFDR > 0.05")) %>%
  ggplot(aes(
    x = LOAD.minus.control, y = -log10(lfdr), 
    color = Significance)) +
  geom_point()+
  theme_minimal() +
  scale_color_manual(values = c("#20854EFF", "grey")) +
  xlab("Effect size") +
  ylab(expression('-log'[10]*'(LFDR)')) +
  xlim(c(-2.5,2.5)) +
  ylim(c(0, 7)) +
  ggtitle("Volcano plot of DMPs")

text_size = 4.5

p.annotated <- p + annotate("segment", x = x_begin, y = yy, xend = x_end, yend = yy, 
           color="black",
           arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
  annotate("text", x = x_begin, size=text_size, hjust=0, y = yy+0.3, label = "Hypermethylated") +
  annotate("segment", x = -x_begin, y = yy, xend = -x_end, yend = yy, 
           color="black",
           arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
  annotate("text", x = -x_end, y = yy+0.3, size=text_size, hjust=0, label = "Hypomethylated") +
  theme(legend.position=c(.2,.75),
        plot.background = element_rect(fill = "white", colour = "white"),
        text = element_text(size=18),
        plot.title = element_text(hjust = 0),
        legend.background = element_rect(color="grey", fill=alpha("white", 1)))



ofile <- file.path(odir, paste0(Sys.Date(), "-volcano-LFDR.png"))
cowplot::save_plot(plot = p.annotated, filename=ofile, 
                   base_width = 7, base_height = 5)


# Pie chart / annotation -----------------------------------------

gr <- df %>% dplyr::filter(lfdr < args$LFDR) %>%
        dplyr::rename(start=chromStart) %>% 
        dplyr::mutate(end=start+1) %>% 
        GRanges()

gene.obj <- readTranscriptFeatures("gencode.hg38.bed")
diff.ann <- annotateWithGeneParts(gr, gene.obj)


# p <- data.frame(diff.ann@annotation) %>% 
#   mutate(Percent = round(diff.ann.annotation)) %>%
#   rownames_to_column("Location") %>%
#   ggplot(aes(fill = Location, values = Percent)) +
#     geom_waffle(n_rows = 10, size = 0.25, color="white") +
#     scale_fill_manual(name = NULL, values = c("#7876B1FF", "#E18727FF", "#FFDC91FF", "#0072B5FF")) +
#   coord_equal() +
#   ggtitle("Location of DMPs") +
#   theme_void() +
#   theme(plot.background = element_rect(fill = "white", colour = "white"),
#         plot.title = element_text(hjust = 0.5, size=20),
#         legend.text=element_text(size=16),
#         legend.position = "bottom",
#         plot.margin = margin(t=5, r=5,b=5,l=5)) +
#   guides(fill=guide_legend(ncol=2))


# ofile <- file.path(odir, paste0(Sys.Date(), "-waffle-DMPs-annotation.png"))
# cowplot::save_plot(plot = p.annotated, filename=ofile, 
#                    base_width = 7, base_height = 5)