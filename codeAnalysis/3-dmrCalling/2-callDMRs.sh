#!/bin/bash

all_files=$(find ../../dataSummaries/dmps-controlLOAD-pc2-ct5 -name *DMPs.forCombP.bed -mmin +5)
# ref_file=../../dataReference/hg38.ncbi.genes.ext.bed


pipeline_wrapper(){
    ifile="${1}"
    prefix="$(dirname $ifile)"
    
    if [ ! -f "$prefix/pvals.acf.bed" ]
    then
        comb-p acf -d 1:1000:100 -c 4 $ifile > "$prefix/acf.txt"
        comb-p slk --acf "$prefix/acf.txt" -c 4 $ifile > "$prefix/pvals.acf.bed"
    else
        print "Skipping $ifile acf step; delete if you want to run"
    fi


    # Some p-values are written as zero because there are no value checks
    # Replace these with 1e-40
    comb-p peaks --dist 500 --seed 0.0001 "$prefix/pvals.acf.bed" \
        |awk '{if ($4 > 0) {print} 
            else print $1"\t"$2"\t"$3"\t1e-40\t"$5
        }' > "$prefix/pvals.regions.bed"

    comb-p region_p -p "$prefix/pvals.acf.bed" \
                -r "$prefix/pvals.regions.bed" \
                -s 100 \
                -c 4 > "$prefix/pvals.regions.tested.bed"
}
    
export -f pipeline_wrapper
parallel --joblog logs/callDMRs.log pipeline_wrapper {} ::: $all_files

# regions_file=../../dataSummaries/DMPs-controlLOAD/combp.regions.bed
# bedtools intersect \
#     -a $regions_file \
#     -b $ref_file > tmp

# prefix="$(dirname $ifile)/combp."

# comb-p pipeline \
#     -c 4 \
#     --seed 1e-5 \
#     --dist 200 \
#     -p $prefix \
#     --region-filter-p 0.01 \
#     $ifile