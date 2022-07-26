library(tidyverse)

i.1 <- "../data/sample-info/DNAMethylation_data_04152022.csv"
i.2 <- "../data/subsamples/2to1_methylation_sample_match_KB.csv"
df.2 <- read_csv(i.2)

df.1 <- read_csv(i.1) %>% separate(sumdx, c("diagnostic_group", "delete")) %>%
        mutate(diagnostic_group = ifelse(diagnostic_group == "Dementia", "LOAD", diagnostic_group)) %>%
        mutate(diagnostic_group = ifelse(diagnostic_group == "Impaired", "Normal", diagnostic_group)) %>%
        mutate(diagnostic_group = ifelse(diagnostic_group == "Normal", "Control", diagnostic_group)) %>%
        rename(sex = gender) %>% 
        mutate(weight = as.numeric(Weight..kg.)) %>%
        mutate(height = as.numeric(Height..cm.)) %>%
        mutate(BMI = 1e4 * weight / height / height)

keep_vars <- c("RA_id", "diagnostic_group", "study_id", "race_primary", "apoe_e1", "apoe_e2", "education",
                "excluded_from_R01", "Charlson.Co.Morbidity.Index.Score", "dxgrp_bin",
                "age_at_visit", "hispanic", "sex", "BMI")

#TODO: Full join with all RA_ids to make master master samplesheet
df <- right_join(df.1, df.2, by = "study_id") %>% dplyr::select(all_of(keep_vars))

set.seed(919)

# LOAD does not get shuffled
ix.load <- which(df$dxgrp_bin == 1)

# Control does
ix <- which(df$dxgrp_bin == 0)
ix <- sample(ix)

ix.1 <- c(ix[1:42], ix.load)
ix.2 <- c(ix[43:84], ix.load)

write_csv(df[ix.1, ], "../data/07-counts-by-chrom-imputed-subset1/master.csv")
write_csv(df[ix.2, ], "../data/07-counts-by-chrom-imputed-subset2/master.csv")