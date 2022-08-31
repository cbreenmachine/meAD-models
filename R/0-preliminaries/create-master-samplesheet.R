# clean-samplesheets.R
# Takes data from three sources
library(tidyverse)
library(readxl)



vars.keep <- c("sample", "cohort", "machine", "pool", 
                "group", "batch", "charlson_score", 
                "bmi", "waist_circumf", "age", "sex", "pack_years",
                "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")

idir <- "../../dataRaw"

design.file <- "sequencingRunInfoFromIllinois/experimental-design.tsv"
mapping.file <- "dataFromArrayAnalysis/fromLigia/sample-IDs-WRAP-ADRC.xlsx"
samples.file <- "dataFromWholeAnalysis/DNAMethylation_data_04152022.csv"

ofile <- "../../dataDerived/masterSamplesheet.csv"

design.df <- read_table(file.path(idir, design.file), show_col_types=F) %>%
    transmute(sample_id = as.character(sample), machine, pool, group, batch)

mapping.df <- read_xlsx(file.path(idir, mapping.file)) %>%
    dplyr::transmute(sample_id = as.character(`RA ID`), source = toupper(source),
        diagnostic_group, source, study_id) %>%
    mutate(study_id = ifelse(!str_detect(study_id, "a"), 
                paste0("adrc", str_pad(study_id, width=5, pad="0")), study_id))

samples.df <- read_csv(file.path(idir, samples.file), show_col_types=F) %>%
        separate(sumdx, c("diagnostic_group", "delete")) %>%
        mutate(diagnostic_group = ifelse(diagnostic_group == "Dementia", "LOAD", diagnostic_group)) %>%
        mutate(diagnostic_group = ifelse(diagnostic_group == "Impaired", "Normal", diagnostic_group)) %>%
        mutate(diagnostic_group = ifelse(diagnostic_group == "Normal", "Control", diagnostic_group)) %>%
        rename(sex = gender) %>% 
        mutate(weight = as.numeric(Weight..kg.),
            height = as.numeric(Height..cm.), 
            BMI = 1e4 * weight / height / height,
            diagnostic_group = toupper(diagnostic_group))


df <- full_join(design.df, mapping.df, by=c("sample_id")) %>%
    full_join(samples.df, by = c("study_id", "source")) %>%
    janitor::clean_names() %>% 
    mutate(diagnostic_group = diagnostic_group_y) %>% 
    dplyr::select(-c("diagnostic_group_y", "diagnostic_group_x")) %>%
    arrange(sample_id)

write_csv(df, ofile)