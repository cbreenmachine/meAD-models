library(tidyverse)

adrc.dir <- "../850k-array"
df <- read_csv("phenos-cleaned-all.csv")
my_files <- list.files(adrc.dir, "adrc") %>% str_remove(".bed")

sub.df <- df %>% filter(adrcnum %in% my_files) %>% select(adrcnum, Alisch)

for (row in 1:nrow(sub.df)){
    old.name <- file.path(adrc.dir, paste0(sub.df[row, 'adrcnum'], ".bed"))
    new.name <- file.path(adrc.dir, paste0(sub.df[row, 'Alisch'], ".bed"))
    print(paste(old.name, new.name))
    file.rename(old.name, new.name)
}
