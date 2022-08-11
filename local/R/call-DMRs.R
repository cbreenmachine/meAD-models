library(data.table)
library(tidyverse)
library(ggbio)
library(wiscR)
# v86 uses hg38
# See https://useast.ensembl.org/info/website/archives/assembly.html
library(EnsDb.Hsapiens.v86)
library(cowplot)
library(ggsci)
library(waffle)
library(patchwork)

#--> Load in data
df <- fread("DMPs.bed", select = c("chrom", "chromStart", "lfdr", "LOAD.minus.control"))
ensdb <- EnsDb.Hsapiens.v86





#--> Top of plot, plots the gene introns/exons
plot_gene_model <- function(ensdb, gene_of_interest){
  
  p <- autoplot(ensdb, ~ symbol == gene_of_interest, 
                stat="reduce", label=F, 
                color="brown", fill="brown") +
    theme_void() +
    ggtitle(paste0(gene_of_interest, " gene model")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab("") +
    scale_x_continuous(position = "top") 
  
  return(p@ggplot)
}




#--> Bottom of plot, scatter of 
plot_points <- function(df, chrom, left, right){
  p <- df %>%
    dplyr::filter(chrom == chrom) %>%
    dplyr::filter(chromStart > left & chromEnd < right) %>%
    dplyr::mutate(Significance = ifelse(lfdr < 0.05, "LFDR < 0.05", "LFDR > 0.05")) %>%
    ggplot(aes(x = chromStart, y = -log10(lfdr) * direction, color = Significance)) +
    geom_point() +
    theme_minimal() +
    xlab("Genomic position") +
    scale_color_nejm() +
    xlim(c(left, right)) +
    ylab(expression('-log'[10]*'(LFDR)')) +
    theme(strip.text.x = element_blank(),
          strip.background = element_rect(colour="white", fill="white"),
          legend.position=c(.9,.87)
    ) 
  
  return(p)
}

# 
# join_via_cow <- function(p1, p2){
#   p <- cowplot::plot_grid(p1, p2, ncol = 1, axis = "b", 
#                           align = "v", rel_heights = c(1, 2)) +
#     theme(plot.background = element_rect(fill = "white", colour = "white"))
#   return(p)
# }


#--> Driver function
make_figure <- function(ensdb, df, gene){
  print(gene)
  subset <- select(ensdb, GeneNameFilter(gene))
  
  left <- min(subset$TXSEQSTART) - pad
  right <- max(subset$TXSEQEND) + pad
  chrom <- paste0("chr", subset$SEQNAME[1])
  
  p1 <- plot_gene_model(ensdb, gene, left, right)
  p2 <- plot_points(df, chrom, left, right)
  p <- join_via_cow(p1, p2)
  
  return(p)
}





#--> PIEZO2 parameters
chrom <- "chr18"
L <- 11.146 * 1e6
R <- 11.149 * 1e6

shade_L <- 11147000
shade_R <- 11148250

p2 <- df %>%
  dplyr::filter(chrom == chrom) %>%
  dplyr::filter(chromStart > L & chromStart < R) %>%
  dplyr::mutate(Type = ifelse(LOAD.minus.control < 0, "Hypermethylated", "Hypomethylated")) %>%
  ggplot(aes(x = chromStart, y = -log10(lfdr), color = Type)) +
  geom_point() +
  theme_minimal() +
  xlab("Genomic position") +
  scale_color_npg() +
  xlim(c(L, R)) +
  ylab(expression('-log'[10]*'(LFDR)')) +
  theme(legend.background = element_rect(color="grey", fill=alpha("white", 1)),
        legend.position=c(.17,.75),
        plot.margin = margin(t=0, r=5,b=5,l=5))
  


#--> All for PIEZO2: Automation coming later


p1 <- plot_gene_model(ensdb, "PIEZO2")

x_range <- layer_scales(p1)$x$range$range
y_below_center <- mean(layer_scales(p1)$y$range$range)


positions <- data.frame(
  x = c(shade_L, x_range[1], x_range[2], shade_R),
  y = c(y_below_center, 0, 0, y_below_center)
)


p1.zoomed <- p1 + 
  geom_polygon(data=positions, aes(x = x, y = y),alpha=0.1) +
  scale_fill_nejm() +
  theme(plot.margin = margin(t=5, r=5,b=0,l=5)) +
  theme(text = element_text(size=16)) +
  theme_void()

ymax <- 2.8

p2.shaded <- p2 +
  annotate("rect", xmin=L, xmax=R,ymin=0,ymax=ymax,alpha=0.05) +
  annotate("rect", xmin=shade_L, xmax=shade_R,
           ymin=0, ymax=ymax, alpha=0.1, color="#E18727FF") +
  annotate("text", x=(shade_L + shade_R) / 2, y=2.5, size=5, color = "#E18727FF",
           label="Differentially Methylated Region") +
  theme(text = element_text(size=16))
# p <- join_via_cow(p1.zoomed, p2.shaded) 
# p
# 
# p <- p1.zoomed / p2.shaded


p <- cowplot::plot_grid(p1.zoomed, NULL, p2.shaded, ncol = 1, axis = "b", 
                        align = "v", rel_heights = c(1, -0.12, 2)) +
  theme(plot.background = element_rect(fill = "white", colour = "white"))

p
cowplot::save_plot(filename="2022-06-17-PIEZO2.png", p, base_width = 8, base_height = 6)








# Volcano Plot ------------------------------------------------------------

ix <- which(sample(df$lfdr) > 0.2)
ix.to.ignore <- ix[1:floor(length(ix) * 0.997)]

x_begin <- 1.2
x_end <- 2.2
yy <- 2.5

p <- df[!ix.to.ignore, ] %>%
  mutate(Significance = ifelse(lfdr < 0.05, "LFDR < 0.05", "LFDR > 0.05")) %>%
  ggplot(aes(
    x = LOAD.minus.control, y = -log10(lfdr), 
    color = Significance)) +
  geom_point()+
  theme_minimal() +
  scale_color_manual(values = c("#20854EFF", "grey")) +
  xlab("Effect size") +
  (expression('-log'[10]*'(LFDR)')) +
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


cowplot::save_plot(plot = p.annotated, filename="2022-06-20-volcano.png", 
                   base_width = 7, base_height = 5)

# Pie chart ---------------------------------------------------------------

loc.df <- data.frame(
  x = rep(1, 4),
  Location = c("Promoter", "Exon", "Intron", "Intergenic"),
  Percent = c(6, 4, 52, 39)  
)

loc.df %>%
  ggplot(aes(fill=Location, x = x, y = Percent)) +
  geom_bar(stat = "identity") +
  scale_fill_nejm() +
  theme_minimal()



percent <- c(6, 4, 52, 38) 
location <- c("Promoter", "Exon", "Intron", "Intergenic")


waffle.plot <- data.frame(percent, location) %>%
  ggplot(aes(fill = location, values = percent)) +
    waffle::geom_waffle(n_rows = 10, size = 0.25, color="white") +
    scale_fill_manual(name = NULL, values = c("#7876B1FF", "#E18727FF", "#FFDC91FF", "#0072B5FF")) +
  coord_equal() +
  ggtitle("Location of DMPs") +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", colour = "white"),
        plot.title = element_text(hjust = 0.5, size=20),
        legend.text=element_text(size=16),
        legend.position = "bottom",
        plot.margin = margin(t=5, r=5,b=5,l=5)) +
  guides(fill=guide_legend(ncol=2))


waffle.plot

cowplot::save_plot(filename = "2022-06-20-waffle.png", waffle.plot)


# png("2022-06-20-DMP-location-waffle.png")
# waffle(z, colors = c("#7876B1FF", "#E18727FF", "#FFDC91FF", "#0072B5FF")) 
# dev.off()