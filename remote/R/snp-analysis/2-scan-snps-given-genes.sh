#!/bin/bash
# 1. Take the 2 input files and put them in bed format,
#  then liftover the second from hg19 to 38

chain=data/hg19ToHg38.over.chain 
ifile1=data/GVC-AD-risk-loci.hg38.csv
ifile2=data/GVC-AD-protective-loci.hg19.csv



if ! [[ -f "${chain}" ]]; then
    wget https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz
    gunzip https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz
    mv hg19ToHg38.over.chain data/
fi


to_bed() {
    ifile="${1}"
    ofile="$(echo ${1} | sed s/csv/bed/)"

    tail -n+3 "${ifile}" \
    | tr ',' ' ' \
    | awk 'BEGIN{OFS="\t"} 
        {print "chr"$2,$3-1,$3}' > "${ofile}"
}

to_bed "${ifile1}"

# Handle the other file
tail -n+3 "${ifile2}"  \
    | sed s/,//4 | sed s/,//4 \
    | sed s/,//4 | sed s/,//4 \
    | sed s/\"// | sed s/\"// \
    | tr "," "\t"| tr ':' '\t' | tr '-' '\t' | tr ' ' '\t' \
    | awk 'BEGIN{OFS="\t"} {print $2,$4,$5,$6}' \
    | cut -f2,3,4 > "$(echo ${ifile2} | sed s/csv/bed/)"
    
../../bin/liftOver \
    data/GVC-AD-protective-loci.hg19.bed \
    data/hg19ToHg38.over.chain \
    data/GVC-AD-protective-loci.hg38.bed \
    unmapped.txt