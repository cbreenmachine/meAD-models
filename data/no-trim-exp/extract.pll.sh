#!/bin/bash
# extract.pll.sh is a parallelized script that converts all bcfs in a directory
# (likely data/03-calls/ into usable bed file s)
# idir : input directory which should contain .bcfs output from gemBS
# odir : where you want to store the output tsv (bed extensions)

idir=.
odir=.

# Extract csv from one bcf (1 and 2 are input and output respectively)
extract() {
    inbcf="$1" # input bcf
    tmpbed="${1}.tmp.tsv"
    outbed="$2" # output
    bcftools query \
        --format '%CHROM\t%POS\t[%CS]\t%REF\t[%MC8{4}]\t[%MC8{5}]\t[%MC8{6}]\t[%MC8{7}]\n' \
        --include 'DP>2 && CG="Y" && DP < 50' "${inbcf}" > "${tmpbed}" 

    python extract_from_bcf.py --ifile "${tmpbed}" --ofile "${outbed}" --merge
}


# List all bcfs in the directory, be mindful of the temp files 
# that are labeled like 107_chr1.bcf using [][][] syntax
files="$(find ${idir} -type f -name "*.bcf")"

# Generate temporary files for parallel
> CRAWLER_INPUT && > CRAWLER_OUTPUT

for infile in ${files}; do
    outfile="${odir}/"$(echo ${infile##*/} | sed -E 's/bcf/bed/')

    if [[ ! -f "${outfile}" ]]
    then
        echo "${infile}" >> CRAWLER_INPUT
        echo "${outfile}" >> CRAWLER_OUTPUT
    fi
done


# GNU parallel needs the function exported
export -f extract 

parallel --link \
    --workdir . \
    --jobs 4 \
    --joblog extract.log \
    extract {1} {2} :::: CRAWLER_INPUT :::: CRAWLER_OUTPUT

# rm CRAWLER_INPUT CRAWLER_OUTPUT
