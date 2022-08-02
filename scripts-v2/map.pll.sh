#!/bin/bash
# TODO: file names for parallel should have sys time appended so we can run multiple instances
idir="${1}"
njobs="${2}"

left_file="map_left_$(basename ${idir})_reads"
right_file="map_right_$(basename ${idir})_reads"

find "${idir}" -name *R1* | sort > ${left_file}
find "${idir}" -name *R2* | sort > ${right_file}

map_wrapper(){
    # Keep in function so it exports
    ref_path=../../reference/GENCODE/gembs-index/h38_no_alt.BS.gem

    # e.g. ../data/pool11/489_R1.fastq.gz
    left="${1}"
    right="${2}"

    # e.g. ../data/pool11/489
    left_root_name=$(echo ${left} | cut -d _ -f1)
    right_root_name=$(echo ${right} | cut -d _ -f1)

    ofile="${left_root_name}.unsorted.sam"
    ofile_alt="{left_root_name}.bam"


    #TODO: reverse these loops; check if json exists

    if [[ $left_root_name == $right_root_name ]]
    then
    echo "Reads match, mapping now"
        if [[ -f "${ofile}" || -f "${ofile_alt}" ]]
        then
            echo "${left} already processed, skipping"
        else
            gem-mapper -I "${ref_path}" \
                --threads 4 \
                --report-file "${left_root_name}.json" \
                --i1 "${left}" --i2 "${right}" -o "${ofile}"
        fi
    else
        echo "Left/Right read names don't match"
    fi
}

log="mapping-$(basename ${idir}).log"

export -f map_wrapper
~/bin/parallel --link --workdir . \
    --jobs "${njobs}" --joblog "${log}" \
    map_wrapper {1} {2} :::: ${left_file} :::: ${right_file}


rm "${left_file}" "${right_file}"
