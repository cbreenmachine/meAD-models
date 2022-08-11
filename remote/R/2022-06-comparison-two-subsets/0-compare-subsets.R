library(data.table)
library(tidyverse)
library(cowplot)
library(viridis)

df.1 <- fread("../../data/07-counts-by-chrom-imputed-subset1/DMPs.bed")
df.2 <- fread("../../data/07-counts-by-chrom-imputed-subset2/DMPs.bed")

alpha <- 0.05

make_locus_col <- function(df){
    df %>% unite("locus", c("chrom", "chromStart"), sep=":") %>% return()
}

df.1.sub <- df.1 %>% make_locus_col
df.2.sub <- df.2 %>% make_locus_col

merged.df <- full_join(df.1.sub, df.2.sub, by = "locus") %>% drop_na

z <- round(cor(merged.df$LOAD.minus.control.x, merged.df$LOAD.minus.control.y), 3)
p <- merged.df %>%
    ggplot(aes(x = LOAD.minus.control.x, y = LOAD.minus.control.y)) +
    # xlim(c(0, 0.9)) +
    # ylim(c(0, 0.9)) +
    geom_hex() +
    scale_fill_viridis() +
    theme_minimal() +
    labs(caption = paste0("Correlation: ", z)) +
    theme(panel.background = element_rect(fill = 'white', colour = 'white'),
        plot.background = element_rect(fill = 'white', colour = 'white'))

cowplot::save_plot(filename="20220730-effectsize-hex.png", p)



z <- round(cor(merged.df$lfdr.x, merged.df$lfdr.y), 3)
p <- merged.df %>%
    filter(lfdr.x < 0.9) %>% filter(lfdr.y < 0.9) %>%
    ggplot(aes(x = lfdr.x, y = lfdr.y)) +
    geom_hex() +
    scale_fill_viridis() +
    theme_minimal() +
    labs(caption = paste0("Correlation: ", z)) +
    theme(panel.background = element_rect(fill = 'white', colour = 'white'),
        plot.background = element_rect(fill = 'white', colour = 'white'))

cowplot::save_plot(filename="20220730-lfdr-hex.png", p)