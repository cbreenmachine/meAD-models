#!/bin/bash

# Constants
chromosomes=$(echo chr{1..22})" chrX chrY"

ifile="${1}"
odir="${2}"

sample=$(basename ${ifile} .bed)

for chr in ${chromosomes};
do
    # For each chromsome, take only those 
    grep -w ${chr} ${ifile} |sed "s/$/\t$sample/" > ${odir}"/"${chr}"."${sample}".bed"
done