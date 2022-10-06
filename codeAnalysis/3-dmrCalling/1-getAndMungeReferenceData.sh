#!/bin/bash

odir=../../dataReference
gff_path="${odir}/hg38.ncbi.ref.gff"
full_bed_path="${odir}/hg38.ncbi.ref.bed"
genes_path="${odir}/hg38.ncbi.genes.bed"
genes_extended_path="${odir}/hg38.ncbi.genes.ext.bed"

# upstream (downstream) nt from gene start (stop)
up=2000 && down=300


#Download data if it does not exist
if [ ! -f $gff_path ]; then
    wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz
    gunzip GRCh38_latest_genomic.gff.gz
    mv GRCh38_latest_genomic.gff $gff_path
fi

# Put everything in a bed for bed tools intersect...
# NOT TESTED YET, CONSIDER ADAPTING BELOW and then GREPING for protein_coding
# if [ -f $full_bed_path ]; then
#     cat $gff_path  \
#         | gff2bed \
#         | sed 's/NC_0000\([1-9][0-9]\)\.[0-9]*\(.*\)/chr\1\2/' \
#         | sed 's/NC_00000\([0-9]\)\.[0-9]*\(.*\)/chr\1\2/' > $full_bed_path
# fi

# Extract just the genes and pack into a bed
if [ ! -f $genes_path ]; then
    # first two lines filter
    # thrid line converts to bed
    # 4and5 clean up gene names a bit
    # last line pads chromosome number to make compatible
    cat $gff_path  \
        | awk '{if ($3 == "gene") {print}}' \
        | grep protein_coding  \
        | gff2bed \
        | sed 's/NC_0000\([1-9][0-9]\)\.[0-9]*\(.*\)/chr\1\2/' \
        | sed 's/NC_00000\([0-9]\)\.[0-9]*\(.*\)/chr\1\2/' > $genes_path
fi


# | sed 's/chr\([1-9]\)[ ]*/chr0\1/' 



# Extent 2kb upstream and 300bp of gene boundaries
# being mindful of strand. Always run this in case
# up/down change
cat $genes_path \
    | awk -v up="$up" -v down="$down" \
        'BEGIN {FS="\t";OFS="\t";}
        {
        if($6 == "+") {print $1,$2-up,$3+down,$4,$5,$6,$7,$8,$9,$10}
        else if($6 == "-") {print $1,$2-down,$3+up,$4,$5,$6,$7,$8,$9,$10}
        }' \
    | awk '$2 >= 0 {print $0}' > $genes_extended_path

# bedtools intersect -a ncbi.genes.2kb300.bed -b DMPs.lfdr05.bed -c > ncbi.geneset.bed
# #-1 in DMPs means hypermethylated
# numTotal=$(wc -l DMPs.lfdr05.bed)
# numHyper=$(bedtools intersect -a ncbi.genes.2kb300.bed -b DMPs.lfdr05.bed -wb | cut -f17 | grep - | wc -l)
