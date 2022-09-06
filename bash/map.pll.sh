#!/bin/bash
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

    # Derived names
    ofile="${left_root_name}.unsorted.sam"
    ojson="${left_root_name}.json"
    obam="${left_root_name}.bam"

    if [ -f "${ofile}" ] || [ -f "${ojson}" ] || [ -f "${obam}" ]; then
        # Output files already exist (json, which reports quality;
        # or bam, which is downstream one step)
        echo "${left} already processed, skipping"
    else
        if [ $left_root_name == $right_root_name ]; then
            
            # Reads are the same name (as they should be)
            echo "Reads match, mapping now"
            
            gem-mapper -I "${ref_path}" \
                --threads 6 \
                --report-file "${left_root_name}.json" \
                --i1 "${left}" --i2 "${right}" -o "${ofile}"
        else
            echo "Left/Right read names don't match"
        fi
    fi
}

log="./logs/mapping-$(basename ${idir}).log"
export -f map_wrapper

if [[ -f "${log}" ]];
then
    echo "Running failed jobs from log file"
    parallel --retry-failed --joblog "${log}" --jobs "${njobs}"
else
    echo "Running from scratch. Won't map files with jsons"
    parallel --link --workdir . \
        --jobs "${njobs}" --joblog "${log}" \
        map_wrapper {1} {2} :::: ${left_file} :::: ${right_file}
fi

rm "${left_file}" "${right_file}"
