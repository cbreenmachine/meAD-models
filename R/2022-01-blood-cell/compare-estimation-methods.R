library(tidyverse)
library(ggsci)
library(cowplot)


vv <- c("sample", "method", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")

df.1 <- read_csv("../../data/study-data/phenos-cleaned.csv") %>%
    mutate(method = "EPIC") %>% select(vv)
df.2 <- read_csv("../../data/07-counts-by-chrom-imputed-subset1//blood-cell-composition.csv") %>%
    mutate(method = "WGMS") %>% select(vv)

df <- rbind(df.1, df.2) %>%
        pivot_longer(cols = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran"),
            values_to = "estimate", names_to = "cell_type") %>%
            pivot_wider(names_from = "method", values_from = "estimate") %>%
            drop_na()

p <- df %>%
        ggplot(aes(x = EPIC, y = WGMS, color = cell_type)) +
        geom_point() +
        theme_minimal() +
        scale_color_nejm() +
        theme(panel.background = element_rect(fill = "white", color = "white"),
            plot.background = element_rect(fill = "white", color = "white")) +
            xlim(c(0,1)) +
            ylim(c(0,1))
cowplot::save_plot(filename = "comparison.png", p)