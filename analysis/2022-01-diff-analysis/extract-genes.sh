#!/bin/bash

ifile=gencode.v40.chr_patch_hapl_scaff.basic.annotation.gtf

awk '{if ($3 == "gene") { print } }' <(tail -n +5  "${ifile}") \
    | grep 'protein_coding' \
    | awk 'BEGIN{OFS="\t"}
        {
            if($7 == "+") {print $1,$4-10000,$5,$14}
            else if($7 == "-") {print $1,$4,$5+10000,$14}
        }' \
    | grep -v 'chrM' | grep 'chr' > TSS10kb.bed

awk '{if ($3 == "gene") { print } }' <(tail -n +5  "${ifile}") \
    | grep 'lncRNA' \
    | awk 'BEGIN{OFS="\t"}{print $1,$4,$5,$14}' \
    | grep -v 'chrM' | grep 'chr' > lncRNA.bed