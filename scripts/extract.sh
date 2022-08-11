#!/bin/bash

# Extract csv from one bcf
ifile="$1"
ofile="$2"

echo ${ifile} ${ofile}

bcftools query \
    --format '%CHROM\t%POS\t[%CS]\t%REF\t[%MC8{4}]\t[%MC8{5}]\t[%MC8{6}]\t[%MC8{7}]\n' \
    --include 'DP>2 && CG="Y" && DP < 50' \
    ${ifile} \
    | python extract_from_bcf.py --ofile ${ofile} --merge