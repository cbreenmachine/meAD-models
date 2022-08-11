#!/bin/bash
idir=../data/00-fastq/

for ifile in $(find "${idir}" -type f -name *fastq.gz);
do
root=$(echo ${ifile} | cut -d "_" -f1)
Rn=$(basename $(echo "${ifile}") | sed -E 's/.*(R[1-2]).*/\1/')
ofile="${root}_${Rn}.fastq.gz"
echo ${ifile} ${ofile}
# mv ${ifile} ${ofile}
done
