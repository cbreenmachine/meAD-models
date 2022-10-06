#!/bin/bash
# Runs like  
# ./5-split-by-chrom.pll.sh ../dataDerived/methylBedBySample/ ../dataDerived/methylBedByChrom/ 20
idir="${1}"
odir="${2}"
njobs="${3}"

chromosomes=$(echo chr{1..22})" chrX chrY"

files="split_$(basename ${idir})_beds"
find "${idir}" -type f -name "*bed" -size +385 -mmin +10 | sort > ${files}

split_wrapper(){
   ifile="${1}"
   odir="${2}"
   chrom="${3}" 

   # e.g. 100
   sample=$(basename ${ifile} .bed)

   # e.g. ../dataOdir/chr22.100.bed
   ofile="${odir}/${chrom}.${sample}.bed"

   if [[ ! -f "${ofile}" ]]; then
      # Adding a sample column, and keep header form original bed
      head -1 ${ifile}| sed "s/$/\tsample/" > "${ofile}"
      grep -w ${chrom} ${ifile} |sed "s/$/\t$sample/" >> "${ofile}"
   fi
}

log="./logs/extract-$(basename ${idir}).log"

export -f split_wrapper

parallel --workdir . \
   --jobs "${njobs}" --joblog "${log}" \
   split_wrapper {1} {2} {3} :::: ${files} ::: ${odir} ::: ${chromosomes}

rm "${files}"
# END
