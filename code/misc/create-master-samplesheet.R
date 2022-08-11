library(tidyverse)
library(readxl)


reid.df <- read_xlsx("../data/study-data/fromLigia/sample-IDs-WRAP-ADRC.xlsx") %>% 
            mutate(RA_id = `RA ID`) %>%
            mutate(source = tolower(source)) 

ix <- which(is.na(reid.df$diagnostic_group))

pad <- function(x){
    y <- paste0("adrc", str_pad(x, 5, "left", "0"))
    return(y)
}

reid.df[ix, "study_id"] <- as.vector(lapply(reid.df[ix, "study_id"], pad))

        
lindsay.df <- read_csv("../data/DNAMethylation_data_04152022.csv")

out.df <- reid.df %>%
    select(RA_id, study_id) %>%
    full_join(lindsay.df, by = "study_id")

write_csv(out.df, "../data/master-samplesheet.csv")


#--> Fully processed files
usable <- tools::file_path_sans_ext(
    tools::file_path_sans_ext(list.files("../data/03-calls", pattern = "bcf.md5")))

tmp.df <- out.df %>%
    filter(RA_id %in% usable) %>%
    select(c(RA_id, study_id, source))

ofile <- "../data/2022-05-16-subset.csv"
write_csv(tmp.df, ofile)