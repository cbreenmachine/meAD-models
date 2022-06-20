#!/bin/bash

ref_file=../data/study-data/experimental-design.tsv

> tmp.txt

for ff in ../data/02-mapping/*bam;
do
    zz=$(basename "${ff##*/}" .bam)
    grep "${zz}" "${ref_file}"  >> tmp.txt
done


cat tmp.txt | tr -s ' ' | cut -d ' '  -f3 | sort | uniq -c

rm tmp.txt
#END

