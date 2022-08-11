#!/bin/bash



map_wrapper() {
    ref_path=../../reference/GENCODE/gembs-index/h38_no_alt.BS.gem
    
    left_file="$1"
    left_root=$(cut -d _ -f1)

    right_file="$2"
    right_root=$(cut -d _ -f1)

    ofile="${left_file%%.*}.sam"
    report_file="${left_file%%.*}.json"

    if [[ "$left_root" == "$right_root" ]]; then
        echo "$ref_path"
        gem-mapper -I "${ref_path}" \
            --threads 4 \
            --report-file "${report_file}" \
            --i1 "${left_file}" --i2 "${right_file}" -o "${ofile}"
    else
        echo "Files do not match"
    fi

}

find 00-fastq/ -name *R1.fastq.gz | sort > LEFT
find 00-fastq/ -name *R2.fastq.gz | sort > RIGHT 

export -f map_wrapper

parallel --jobs 2 --link map_wrapper :::: LEFT :::: RIGHT