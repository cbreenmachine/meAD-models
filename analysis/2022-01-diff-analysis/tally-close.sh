#!/bin/bash
ifile=chr18/smooth-150-PCs-2/DMPs.bed
ofile=chr18-TSS10kb.bed

bedtools intersect -a <(tail -n +2 "${ifile}") \
    -b TSS10kb.bed > "${ofile}"

wc -l ${ifile}
uniq ${ofile} | wc -l 


bedtools intersect -a <(tail -n +2 "${ifile}") \
    -b lncRNA.bed > "chr18-lncRNA.bed"
