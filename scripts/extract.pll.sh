#!/bin/bash
# extract.pll.sh is a parallelized script that converts all bcfs in a directory
# (likely data/03-calls/ into usable bed file s)
# idir : input directory which should contain .bcfs output from gemBS
# odir : 

# prefix=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
prefix="$(dirname $0)/" 
echo ${prefix}

idir="${1}"
odir="${2}"

# Extract csv from one bcf
extract() {
    bcftools query \
        --format '%CHROM\t%POS\t[%CS]\t%REF\t[%MC8{4}]\t[%MC8{5}]\t[%MC8{6}]\t[%MC8{7}]\n' \
        --include 'DP>2 && CG="Y" && DP < 50' "$1" \
        |python "${prefix}extract_from_bcf.py" --ofile "$2" --merge
}

# List all bcfs in the directory, be mindful of the temp files 
# that are labeled like 107_chr1.bcf using [][][] syntax
files="$(find ${idir} -type f -name "*[0-9][0-9][0-9].bcf")"

# Generate temporary files for parallel
> CRAWLER_INPUT && > CRAWLER_OUTPUT

for infile in ${files}; do
    outfile=$(echo ${infile} | sed -E 's/bcf/bed/')
    echo "${infile}" >> CRAWLER_INPUT
    echo "${outfile}" >> CRAWLER_OUTPUT
done

export -f extract 

parallel --link  \
    --workdir . \
    --S nebula-2,nebula-3,nebula-4 \
    --tempdir ./ --jobs 50% \
    --joblog extract.log \
    extract {1} {2} :::: CRAWLER_INPUT :::: CRAWLER_OUTPUT

# rm CRAWLER_INPUT CRAWLER_OUTPUT