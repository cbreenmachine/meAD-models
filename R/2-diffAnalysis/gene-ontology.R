library(clusterProfiler)
library(tidyverse)
library(data.table)

pull_gene_id <- function(s){
    str_extract(s, "GeneID:[0-9]+") %>%
    str_remove("GeneID:") %>%
    return()
}

pull_gene_id_v <- Vectorize(pull_gene_id)

df <- fread("ncbi.geneset.bed") %>%
    dplyr::mutate(gene_id = pull_gene_id_v(V10))


df.2 <- df %>% transmute(gene_id, count = V11) %>% filter(count >= 5) 
gene_list <- as.vector(df.2$gene_id)

go.out <- clusterProfiler::enrichGO(gene_list, OrgDb = 'org.Hs.eg.db', ont="ALL", pAdjustMethod="fdr" )
go.out

go.out %>% group_by(ONTOLOGY) %>% summarize(N = n())

tmp <- dplyr::filter(go.out, ONTOLOGY == "BP")
tmp$Description


DMPs <- fread("DMPs.lfdr05.bed") %>% mutate(direction = ifelse(direction == -1, "Hyper", "Hypo"))

N <- nrow(DMPs)

DMPs %>% group_by(direction) %>%
    summarize(Perc = n() / N) 