#!/bin/bash
# call.pll.sh
# Operates on directory of (usually 20) files and takes input bam files and outputs
# bcf variant call files
# Runs like `./call.pll.sh ../data/pool17/ 10`
# Runs in gem environment--need parallel, bcftools, and samtools
# Indexing takes ~10-15 mins
idir="${1}"
njobs="${2}"

files="extract_$(basename ${idir})_bcfs"
find "${idir}" -type f -name "*bcf" -size +14G -mmin +30 | sort > ${files}

extract_wrapper(){
    # e.g. ../data/pool01/110.bcf
    ifile="${1}"
    odir=../dataDerived/methylBedBySample/

    # Derive output name by replacing bcf with bed, 
    # swapping out directory name
    ofile=$(echo "${odir}$(basename $ifile bcf)bed")
    # tmpfile=$(echo "${odir}$(basename $ifile bcf)tmp.tsv")
   
    if ! [[ -f "${ofile}" ]]; then
        bcftools query \
            --format '%CHROM\t%POS\t[%CS]\t%REF\t[%MC8{4}]\t[%MC8{5}]\t[%MC8{6}]\t[%MC8{7}]\n' \
            --exclude 'DP<=2 | DP>=50 | CG="N" | CG="H" | CG="?"' "${ifile}" \
        | python extract_from_bcf.py --ofile "${ofile}" --merge
    else
        echo "${ofile} already exists, delete if you want to process"        
    fi  
}

log="./logs/extract-$(basename ${idir}).log"

export -f extract_wrapper
parallel --link --workdir . \
    --jobs "${njobs}" --joblog "${log}" \
    extract_wrapper :::: ${files}

rm "${files}"