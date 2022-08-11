#!/bin/bash

idir="../../data/pool01"
regions_file="data/GVC-AD-protective-loci.hg38.bed"

bcfs_file=BCFS

find "${idir}" -name "*bcf" > "${bcfs_file}"


filter_wrapper() {
    ifile="${1}"
    # ofile="${}"

    bcftools query \
        --include 'FILTER="PASS"' \
        --regions "${regions_file}" \
        --format '%CHROM\t%POS\t%REF\t%ALT\n' \
        "${1}"  # | grep -v "0/0:"

}

export -f filter_wrapper
parallel --jobs 2 \
    --workdir . \
    filter_wrapper {1} :::: BCFS