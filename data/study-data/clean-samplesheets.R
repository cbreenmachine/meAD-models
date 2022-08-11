library(tidyverse)

tmp <- read_tsv("ADRC-samplesheet-cleaned.tsv")
pheno.df <- read_csv("archived/master-WB.csv") %>% 
            inner_join(tmp, by = c("adrcnum" = "ADRC_ID")) %>% 
            mutate(sample = as.character(ALISCH_ID))
design.df <- read_table("experimental-design.tsv") %>% 
            dplyr::mutate(sample = as.character(sample))

df <- inner_join(pheno.df, design.df, by = "sample") 

write_csv(df, "phenos-cleaned-all.csv")

df.filt <- df %>%
    rename(sample = Alisch, batch = batch.y) %>%
    dplyr::select(c(sample, cohort, machine, pool, group, batch,
                    CD8T, CD4T, NK, Bcell, Mono, Gran,
                    charlson_score, bmi, waist_circumf, 
                    age, sex, pack_years))

write_csv(df.filt, "phenos-cleaned.csv")