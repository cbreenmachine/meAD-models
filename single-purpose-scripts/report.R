library(tidyverse)

df <- read_csv("../data/master-samplesheet.csv")

include <- tools::file_path_sans_ext(list.files("../data/02-mapping/"))
include <- tools::file_path_sans_ext(list.files("../data/03-calls/", pattern="[0-9][0-9][0-9].bcf$"))

tmp <- df %>% filter(RA_id %in% include) %>% select(study_id, RA_id, source)

write_csv(tmp, "2022-05-10-subset-analysis.csv")

df %>% filter(RA_id %in% include) %>%
    group_by(diagnostic_group) %>%
    summarize(count = n())
