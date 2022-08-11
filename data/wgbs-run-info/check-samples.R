library(tidyverse)
library(readxl)

ofile <- "summary.txt"

uw.files <- list.files(pattern = "ligia")
il.files <- list.files(pattern = "Report")

uw.df <- bind_rows(lapply(uw.files, function(x) select(readxl::read_xlsx(x, col_types = "text"), c(batch, sample))))




sheetNr <- lapply(il.files, function(x) length(excel_sheets(x)))

lapply(1:length(il.files), function(x) {
    read_xlsx(il.files[x], sheet = sheetNr[[x]] ) %>%
    drop_na(Sample) %>%
    select(c(Project, Sample))
}) %>% bind_rows() %>% 
    mutate(Project = str_to_title(word(Project, sep = "_", 2))) -> il.df


while(length(ix <- which(is.na(il.df$Project))) > 0){
  il.df$Project[ix] <- il.df$Project[ix -1]
}



tmp <- paste("Number of unique sample identifiers from Illinois:", length(unique(il.df$Sample)))
write.table(tmp, file = ofile, col.names = F, row.names = F, quote = F)


tmp <- paste("Sample in UW samplesheets but not in Illinois Data:", setdiff(uw.df$sample, il.df$Sample)[-3])
write.table(tmp, file = ofile, append = TRUE,  quote = FALSE, col.names = F, row.names = F)
tmp <- paste("Sample in Illinois data but not in UW samplesheets:", setdiff(il.df$Sample, uw.df$sample))
write.table(tmp, file = ofile, append = TRUE,  quote = FALSE, col.names = F, row.names = F)


write.table("Size of each pool:", file = ofile, append = TRUE,  quote = FALSE, col.names = F, row.names = F)
il.df %>%
    group_by(Project) %>%
    summarize(Count = n()) %>%
    write.table(file = ofile, append = TRUE,  quote = FALSE,col.names = F, row.names = F)


write.table("Location of 20Xs:", file = ofile, append = TRUE,  quote = FALSE,col.names = F, row.names = F)
il.df %>% 
    filter(str_detect(Sample, "^20")) %>%
    write.table(file = ofile, append = TRUE,  quote = FALSE,col.names = F, row.names = F)


