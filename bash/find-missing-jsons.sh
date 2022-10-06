#! /bin/bash

find ../data/ -name *am |grep -v tmp | cut -d "/" -f4 | cut -d "." -f1 |sort|uniq | wc -l

for ff in $(find ../data -name *unsorted*)
do 
    json=$(echo $ff | sed s/unsorted\.sam/json/)

    if [ ! -f ${json} ]; then echo $ff; fi

done