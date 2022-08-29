#!/bin/bash
# bash generate-DMP-report.sh DMPs.lfdr01.bed &> report.lfdr01.txt

reffile=genes2kb300.bed
ifile="${1}"
ofile="report.${1}"

echo "Total DMPs : $(wc -l ${ifile})" > ${ofile}
echo "Total DMPs within 2kb of TSS / 300bp of end : $(bedtools intersect -a ${reffile} -b ${ifile} | wc -l)" >> ${ofile}

tmp=$(bedtools intersect -a ${ifile} -b ${reffile} -c | awk '$8 > 1' | wc -l)
echo "Number of DMPs overlapping multiple upstream regions : ${tmp}" >> ${ofile}

# Print list of genes withs 5 or more DMPs per 
echo "---------------------------------" >> ${ofile}
echo "Genes with 5 or more DMPs : " >> ${ofile}
bedtools intersect -a ${reffile} -b ${ifile} -c | awk '$6 >= 3' >> ${ofile}

bedtools intersect -a ${reffile} -b <(cut ${ifile} -f1-3,7) -wb | cut -f5,9 > tmp.bed


echo "---------------------------------" >> ${ofile}
echo "Direction of effect per DMR : " >> ${ofile}

> signs.bed
for gene in $(cut tmp.bed -f1 | uniq)
do
tmp=$(grep ${gene} tmp.bed | cut -f2)
numTotal=$(echo ${tmp} | tr -cd '1' | wc -c)
numNeg=$(echo ${tmp} | tr -cd '-' | wc -c)
echo -e ${gene}" \t"${numTotal}" \t"${numNeg} >> signs.bed
done

echo -e "Gene\tTotalDMPs\tNumPositive" >> ${ofile}
awk '$2 >= 3$' signs.bed >> ${ofile}