library(tidyverse)

df <- read_csv('../../data/sample-info/DNAMethylation_data_04152022.csv') %>% 
        separate(sumdx, c("diagnostic_group", "delete")) %>%
        mutate(diagnostic_group = ifelse(diagnostic_group == "Dementia", "LOAD", diagnostic_group)) %>%
        mutate(diagnostic_group = ifelse(diagnostic_group == "Impaired", "Normal", diagnostic_group)) %>%
        mutate(diagnostic_group = ifelse(diagnostic_group == "Normal", "Control", diagnostic_group)) 

df %>%
group_by(diagnostic_group, gender) %>%
summarize(Count = n()) %>%
write_csv("sex-diagnosis.csv")