#!/bin/bash

dir="$1"

echo "Working on $1"

#TODO: clean filename function

sort_wrapper(){
    ifile="${1}"
    ofile=$(echo ${ifile} | sed 's/unsorted\.sam/bam/')
    echo ${ifile} ${ofile}

    # samtools sort -@4 -m4g "${ifile}" -o "${ofile}"
}



