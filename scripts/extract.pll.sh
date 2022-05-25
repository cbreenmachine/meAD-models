#!/bin/bash
# extract.pll.sh is a parallelized script that converts all bcfs in a directory
# (likely data/03-calls/ into usable bed file s)
# idir : input directory which should contain .bcfs output from gemBS
# odir : 

idir="${1}"
odir="${2}"

# mkdir -p "${odir}"

# Extract csv from one bcf
extract() {
    ifile="$1"
    ofile="$2"

    bcftools query \
        --format '%CHROM\t%POS\t[%CS]\t%REF\t[%MC8{4}]\t[%MC8{5}]\t[%MC8{6}]\t[%MC8{7}]\n' \
        --include 'DP>2 && CG="Y" && DP < 50' \
        ${ifile} \
        | python extract_from_bcf.py --ofile ${ofile} --merge
}

# List all bcfs in the directory, be mindful of the temp files 
# that are labeled like 107_chr1.bcf
files="$(find ${idir} -type f -name "*[0-9][0-9][0-9].bcf")"

# Generate temporary files for parallel
> CRAWLER_INPUT && > CRAWLER_OUTPUT

for infile in ${files}; do
    outfile=$(echo ${infile} | sed -E 's/03-calls/04-extract/' | sed -E 's/bcf/bed/')
    echo "${infile}" >> CRAWLER_INPUT
    echo "${outfile}" >> CRAWLER_OUTPUT
    
done

export -f extract 
#--workdir . 
#TODO:Change workdir so that 
parallel --link \
    --tempdir ./ --jobs 8 \
    --joblog extract.log \
    extract {1} {2} :::: CRAWLER_INPUT :::: CRAWLER_OUTPUT

rm CRAWLER_INPUT CRAWLER_OUTPUT