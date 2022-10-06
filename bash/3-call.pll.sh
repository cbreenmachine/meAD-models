#!/bin/bash
# call.pll.sh
# Operates on directory of (usually 20) files and takes input bam files and outputs
# bcf variant call files
# Runs like `./call.pll.sh ../data/pool17/ 10`
# Runs in gem environment--need parallel, bcftools, and samtools
# Indexing takes ~10-15 mins

idir="${1}"
njobs="${2}"

export PATH="~/miniconda3/envs/gem/bin/:$PATH"

files="call_$(basename ${idir})_bams"
find "${idir}" -type f -name [0-9][0-9][0-9].bam | sort > ${files}

call_wrapper(){

    # Keep in function so it exports
    ref_path=../../reference/GENCODE/h38_no_alt.fa

    # e.g. ../data/pool11/489_R1.fastq.gz
    ifile="${1}"
    ofile=$(echo ${ifile} | sed s/bam/bcf/)
   
    if ! [[ -f "${ofile}" ]]; then
        # Output DNE --> run calling
        # -p (for paired end) flag is used in github readme.
        #but is not a valid option when actually running bs_call
        samtools view -h -f 2 "${ifile}" \
            | bs_call --reference "${ref_path}" -L4 \
            | bcftools view -Ob - > "${ofile}"

        # Also index the output bcf
        bcftools index "${ofile}"
    else
        echo "${ofile} already exists, delete if you want to process"        
    fi
}

log="./logs/calling-$(basename ${idir}).log"

export -f call_wrapper
parallel --link --workdir . \
    --jobs "${njobs}" --joblog "${log}" \
    call_wrapper {1} :::: ${files}
    
rm "${files}"
#END