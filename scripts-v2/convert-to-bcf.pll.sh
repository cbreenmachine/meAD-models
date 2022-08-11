#!/bin/bash
idir="${1}"
njobs="${2}"

# Where we store the files to be sorted
sort_file="convert_$(basename ${idir})"
find "${idir}" -type f -name "*vcf" | sort > "${sort_file}"


convert_wrapper() {
    ifile="${1}"
    ofile=$(echo ${ifile}| sed s/vcf/bcf/)

    echo ${ifile}
    echo ${ofile}
    bcftools view -Ob "${ifile}" > "${ofile}"
    bcftools index "${ofile}"
}


export -f convert_wrapper
~/bin/parallel --link --workdir . \
    --jobs "${njobs}" --joblog "${log}" \
    convert_wrapper :::: "${sort_file}"

rm "${sort_file}"