#!/bin/bash
# prepare-entrez.sh
# Need gene IDs from NCBI to do anlaysis in clusterProfiler in R

# wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz
# gunzip GRCh38_latest_genomic.gff.gz

cat GRCh38_latest_genomic.gff  \
    | awk '{if ($3 == "gene") {print}}' \
    | grep protein_coding  \
    | gff2bed \
    | sed 's/NC_0000\([1-9][0-9]\)\.[0-9]*\(.*\)/chr\1\2/' \
    | sed 's/NC_00000\([0-9]\)\.[0-9]*\(.*\)/chr\1\2/' > ncbi.genes.bed

wc -l ncbi.genes.bed

awk 'BEGIN {FS="\t";OFS="\t";}
    {
        if($6 == "+") {print $1,$2-2000,$3+300,$4,$5,$6,$7,$8,$9,$10}
        else if($6 == "-") {print $1,$2-300,$3+2000,$4,$5,$6,$7,$8,$9,$10}
    }' ncbi.genes.bed > ncbi.genes.2kb300.bed

# Grep is hacky, seehttps://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_assembly_report.txt
# To map random contigs
   

bedtools intersect -a ncbi.genes.2kb300.bed -b DMPs.lfdr05.bed -c > ncbi.geneset.bed

#-1 in DMPs means hypermethylated
numTotal=$(wc -l DMPs.lfdr05.bed)
numHyper=$(bedtools intersect -a ncbi.genes.2kb300.bed -b DMPs.lfdr05.bed -wb | cut -f17 | grep - | wc -l)