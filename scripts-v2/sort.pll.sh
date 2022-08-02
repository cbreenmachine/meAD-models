#!/bin/bash

idir="${1}"
njobs="${2}"

# Where we store the files to be sorted
sort_file="sort_$(basename ${idir})"
find "${idir}" -type f -name "*unsorted.sam" | sort > "${sort_file}"


#TODO: clean filename function

sort_wrapper(){
    ifile="${1}"
    ofile=$(echo ${ifile} | sed 's/unsorted\.sam/bam/')
    json=$(echo ${ifile} | sed 's/unsorted\.sam/json/')
    

    # First check if ofile exists, if it does no need to sort
    if ! [[ -f "${ofile}" ]]; then
        
        # Then check if json is there (this means mapping ran successfully)
        if [[ -f "${json}" ]]; then
            echo ${ifile} ${ofile}
            samtools sort -@4 -m4g "${ifile}" -o "${ofile}"
        else
            echo "${ofile} does not have json; re-map if you want to run sorting"
        fi
    else
        echo "${ofile} already exists; delete if you want to run sorting"
    fi
    
    

    
}


log="sorting-$(basename ${idir}).log"

export -f sort_wrapper

parallel --link --workdir . \
    --jobs "${njobs}" --joblog "${log}" \
    sort_wrapper {1} :::: "${sort_file}"

rm "${sort_file}"
