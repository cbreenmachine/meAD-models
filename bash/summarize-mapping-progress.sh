#!/bin/bash

ref_file=../dataDerived/masterSamplesheet.csv

tmp_ref_file=./tmpRef.csv

cat ${ref_file} | cut -d ',' -f1,3 | sort | uniq > ${tmp_ref_file}

> processed.txt
for zz in $(find ../data/ -name *bam | rev | cut -d '/' -f1 | rev | cut -d '.' -f1 | sort |uniq);
do
    grep "${zz}" "${tmp_ref_file}"  >> processed.txt
done

cat processed.txt | cut -d ','  -f2 | sort | uniq -c

#> processed.txt
#for zz in $(find ../data/ -name *am | rev | cut -d '/' -f1 | rev | cut -d '.' -f1 | sort |uniq);
#do
#    grep "^${zz}," "${ref_file}"  >> processed.txt
#done
#cat processed.txt  |cut -d  ","  -f26 | sort | uniq -c

#rm processed.txt

#END
