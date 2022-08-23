#!/bin/bash

ref_file=../data/study-data/experimental-design.tsv

> processed.txt
for zz in $(find ../data/ -name *bcf | rev | cut -d '/' -f1 | rev | cut -d '.' -f1 | sort |uniq);
do
    grep "${zz}" "${ref_file}"  >> processed.txt
done

cat processed.txt | tr -s ' ' | cut -d ' '  -f3 | sort | uniq -c

ref_file=../data/sample-info/master-samplesheet.csv

> processed.txt
for zz in $(find ../data/ -name *bcf | rev | cut -d '/' -f1 | rev | cut -d '.' -f1 | sort |uniq);
do
    grep "^${zz}," "${ref_file}"  >> processed.txt
done
cat processed.txt  |cut -d  ","  -f26 | sort | uniq -c

rm processed.txt

#END
