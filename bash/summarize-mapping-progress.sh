#!/bin/bash

ref_file=../dataDerived/masterSamplesheet.csv

> processed.txt
for zz in $(find ../data/ -name *am | rev | cut -d '/' -f1 | rev | cut -d '.' -f1 | sort |uniq);
do
    grep "${zz}" "${ref_file}"  >> processed.txt
done

cat processed.txt | tr -s ' ' | cut -d ' '  -f3 | sort | uniq -c

> processed.txt
for zz in $(find ../data/ -name *am | rev | cut -d '/' -f1 | rev | cut -d '.' -f1 | sort |uniq);
do
    grep "^${zz}," "${ref_file}"  >> processed.txt
done
cat processed.txt  |cut -d  ","  -f26 | sort | uniq -c

rm processed.txt

#END
