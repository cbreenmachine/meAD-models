#!/bin/bash
ifile=data/gencode.v40.chr_patch_hapl_scaff.basic.annotation.gtf
ofile=regions.tsv

# GTFs are one-based, and regions of bctools defaults to that, so no problem

if ! [[ -f "${ofile}" ]]; then
awk '{if ($3 == "gene") { print } }' <(tail -n +5  "${ifile}") \
    | grep 'protein_coding' \
    | grep '"PSEN1"\|"PSEN2"\|"APP"\|"APOE"\|"RYR1"\|"RYR2"\|"RYR3"' \
    | cut -f1,4,5 > "${ofile}"
fi

filter_wrapper() {
    ifile="${1}"
    # ofile="${}"

    bcftools view \
        --include 'FILTER="PASS" & QUAL > 50' \
        --regions-file regions.tsv \
        -m2 -M2 -v snps "${ifile}" \
    | bcftools query \
        --format '%CHROM\t%POS\t%REF\t%ALT\t%QUAL[\t%GT]\n' 
}

filter_wrapper ../../data/pool01/101.bcf 