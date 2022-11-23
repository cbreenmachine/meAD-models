library(tidyverse)

ref.genes <- readxl::read_xlsx("../../dataReference/2022-natureGenetics.xlsx") %>% pull(query)
dm.genes <- read.table("../../figs/2022-paper/lfdr0.05-effect0/genesTestedSig.txt")$V1


intersect(ref.genes, dm.genes)
