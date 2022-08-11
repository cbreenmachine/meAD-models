# cobble-reports.R 
# Crawls through data directory, looking for QC reports output by gemBS, and writes out two data frames;
library(rvest)
library(tidyverse)
library(janitor)

ofile <- "quality.csv"
array_samples <- as.character(pull(read_csv("../../data/meta/array-samples.csv"), "sample"))

#--> Pull in viable directories on which to operate
all_files <- list.files(path = "../../data/", pattern = "index.html", full.names = TRUE, recursive=TRUE)
call_files <- all_files[str_detect(all_files, "calling")]
map_files <- all_files[str_detect(all_files, "mapping")]

#--> Functions
cobble_mapping <- function(path){
  path %>% 
    read_html() %>% 
    html_nodes(xpath = '//*[@id="hor-zebra"][1]') %>%
    .[[1]] %>% 
    html_table(fill = TRUE) %>% 
    mutate(Sample = str_remove(Sample, "» \n"),
           Unique = as.numeric(word(`Unique (%)`))) %>% 
    clean_names() %>%
    mutate(unique_percent = gsub("[()]", "", word(unique_percent, 2))) %>%
    select(-c(conv_rate, over_conv_rate)) %>%
    return()
}


cobble_calling <- function(path){
    # Takes a calling/index.html and outputs dataframe with 
  path %>% 
    read_html() %>% 
    html_nodes(xpath = '//*[@id="hor-zebra"][1]') %>%
    .[[1]] %>% 
    html_table(fill = TRUE) %>% 
    clean_names() %>%
    mutate(sample = as.character(sample),
        aligned = word(aligned), 
        uniquely_aligned = word(uniquely_aligned, 1),
        passed = gsub("[()]", "", word(passed, 1)),
        passed_variants = gsub("[()]", "", word(passed_variants, 2)),
        median_cpg_cov = str_remove(med_cp_g_cov, "x"),
        med_cov_passed_variants = str_remove(med_cov_passed_variants, "x")) %>%
    select(-c(reports, med_cp_g_cov)) %>%
    return()
}

#--> Pull it all together
map.df <- do.call(rbind, lapply(map_files, cobble_mapping))
call.df <- do.call(rbind, lapply(call_files, cobble_calling))
df <- inner_join(map.df, call.df, by = "sample") %>% 
        mutate(gigabases = 150 * reads / 1e9) %>%
        filter(sample %in% array_samples)


write_csv(df, ofile)