#!/bin/bash

idir="${1}"

find "${idir}" -name *R1* | sort > map_left_reads
find "${idir}" -name *R2* | sort > map_right_reads

map_wrapper(){
    # Keep in function so it exports
    ref_path=../../reference/GENCODE/gembs-index/h38_no_alt.BS.gem

    # e.g. ../data/pool11/489_R1.fastq.gz
    left="${1}"
    right="${2}"

    # e.g. ../data/pool11/489
    left_root_name=$(echo ${left} | cut -d _ -f1)
    right_root_name=$(echo ${right} | cut -d _ -f1)
   
    if [[ $left_root_name == $right_root_name ]]
    then
    echo "Reads match, mapping now"
    gem-mapper -I "${ref_path}" \
        --threads 4 \
        --report-file "${left_root_name}.json" \
        --i1 "${left}" --i2 "${right}" -o "${left_root_name}.unsorted.sam"
    else
    echo "Left/Right read names don't match"
    fi
}

export -f map_wrapper
parallel --link --workdir . \
    --jobs 5 \
    map_wrapper {1} {2} :::: map_left_reads :::: map_right_reads
