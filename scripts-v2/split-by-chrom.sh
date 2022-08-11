#!/bin/bash

# Constants
chromosomes=$(echo chr{1..22})" chrX chrY"

ifile="${1}"
odir="${2}"

sample=$(basename ${ifile} .bed)


for chr in ${chromosomes};
do
    # For each chromsome, take only those 
    ofile=${odir}"/"${chr}"."${sample}".bed"

    # Adding a sample column, and keep header form original bed
    head -1 ${ifile}| sed "s/$/\tsample/" > "${ofile}"
    grep -w ${chr} ${ifile} |sed "s/$/\t$sample/" >> "${ofile}"

done

echo "${ifile} complete"

# END
