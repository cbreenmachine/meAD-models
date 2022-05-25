#!/bin/bash

source configuration.txt

prepare() {
  cd ${idir}
  gemBS prepare -c ${CONF_OUT} -t ${META_OUT}
}

map() {
  # MAPPING
  cd ${idir}
  parallel -S nebula-5,nebula-2,mastodon-3 --joblog ${date}-map.log --nonall --workdir . gemBS map --no-merge
  cd ${home}
  # END MAPPING
}

call(){
  # Calling
  cd ${idir}
  parallel -S ${RUN_SERVERS} --joblog ${date}-call.log --nonall --workdir . gemBS call
  gemBS report
  cd ${home}
  # END CALLING
}



extract() {
  # Extract csv from one bcf
  ifile="$1"
  ofile="$2"

  bcftools query --format '%CHROM\t%POS\t[%CS]\t%REF\t[%MC8{4}]\t[%MC8{5}]\t[%MC8{6}]\t[%MC8{7}]\n' \
    --include 'DP>2 && CG="Y" && DP < 50' ${ifile} | python extract_from_bcf.py --ofile ${ofile} --merge
}

clean_trim_dir() {

  cd ${idir}
  cd ${FASTQ_TRIMMED_PATH}
  # remove anything that didn't finish
  rm *trimmed.fq

  unique_nums=$(find . -name "*.txt" | cut -c3-5) # get list like 235 236 237 ...
  echo ${unique_nums}
  for un in ${unique_nums}; 
  do 
    if [[ $(ls ${un}* | wc -l) -lt 4 ]]; then
        echo "Deleting files starting with ${un}"
        rm ${un}*
    fi
  done
  cd ${home}

}


trim() {
    
    > LEFT # LEFT reads stored line-by-line in this file
    > RIGHT
    for ff in $(find ${idir}${FASTQ_PATH} -name *R1.fastq.gz); do 
        un=$(basename ${ff} | cut -c1-3) # grab the 123 from 00-fastq/123_stuff.fastq.gz

        # a successful trim_galore run produces 4 files starting with 123
        # if less than 4, needs trimming
        if [[ $(find "${FASTQ_TRIMMED_PATH}" -name "${un}*" | wc -l) -lt 4 ]] 
        then
            echo "${ff}" >> LEFT; 
            echo ${ff} | sed -e s/R1/R2/ | sed -e s/val_1/val_2/ >> RIGHT; 
        fi
    done

    cat LEFT 
    cat RIGHT

    # --link creates a mapping between the lines in LEFT and lines in RIGHT 
    # (one-to-one like python's zip instead of pairwise combinations)
    # the fourth ':' means cat LEFT and RIGHT (don't treat as variable/expansion)
    parallel --link -S nebula-5 \
      --jobs 3 --workdir . --joblog trim.log \
        trim_galore --phred33 --cores 6 \
        --output_dir ${FASTQ_TRIMMED_PATH} \
        --dont_gzip --paired {1} {2} :::: LEFT :::: RIGHT
    rm LEFT RIGHT # cleanup

    cd ${home}
}

extract_all() {

  files="$(find ${data_dir} -type f -name "*[0-9][0-9][0-9].bcf")"
  > CRAWLER_INPUT && > CRAWLER_OUTPUT

  for infile in ${files}; do
      outfile=$(echo ${infile} | sed -E 's/03-calls/04-extract/' | sed -E 's/bcf/bed/')
      echo "${infile}" >> CRAWLER_INPUT
      echo "${outfile}" >> CRAWLER_OUTPUT
      new_dir=$(dirname "${outfile}")
      mkdir -vp "${new_dir}"
  done

  echo ""
  echo "Will extract methylation on the following files:"
  cat CRAWLER_INPUT

  export -f extract
  parallel --link --workdir . --tempdir ./ --jobs 8 --joblog extract.log \
      extract {1} {2} :::: CRAWLER_INPUT :::: CRAWLER_OUTPUT

  rm CRAWLER_INPUT CRAWLER_OUTPUT
}


regions=$(echo chr{1..22})" chrX chrY"

split_by_chr() {
  idir="${1}"
  odir="${2}"

  # First create the files that will be called as 
  # odir/chr1.tsv (not they're no longer realy BED files)
  for rr in ${regions};
  do
    > "${odir}/${rr}.tsv"
  done


  for ff in $(find ${idir} -name *bed);
  do
    for rr in ${regions};
    do
      grep -w $rr <(tail -n +2  ${ff}) >> "${odir}/${rr}.tsv"
    done
  done

}


# if invoked as a script rather than sourced, call function named on argv via the below;
# note that this must be the first operation other than a function definition
# for $_ to successfully distinguish between sourcing and invocation:
#[[ $_ != $0 ]] && return

# make sure we actually *did* get passed a valid function name
if declare -f "$1" >/dev/null 2>&1; then
  # invoke that function, passing arguments through
  "$@" # same as "$1" "$2" "$3" ... for full argument list
else
  echo "Function $1 not recognized" >&2
  exit 1
fi







# prepare() {
#   cd ${idir}
#   # BEGIN MAKE META FILE
#   echo "barcode,dataset,end1,end2" > ${META_OUT}
#   for f in "${FASTQ_TRIMMED_PATH}"*val_1.fq
#   do
#       # f="../../data/something.csv"
#       # --> ${f##*/} is something.csv
#       barcode=$(echo ${f##*/} | cut -d "_" -f 1 | sed -E 's/R/s/')
#     # dataset=$(echo ${f##*/} | sed -E 's/val_[1-2].fq//' | sed -E 's/_R1_//')
#       file1=$(echo ${f##*/})
#       file2=$(echo ${f##*/} | sed -E 's/R1/R2/' | sed -E 's/val_1/val_2/')
#       echo "${barcode},${barcode},${file1},${file2}" >> ${META_OUT}
#   done
#   # END MAKE META FILE
#   echo "Preparing gemBS, should only take a second"
#   echo "${CONF_TXT}" > ${CONF_OUT}
#   # GEMBS PREPARATION AND CONSOLE OUTPUT
#   rm -rf .gemBS # clears some hanging errors
#   gemBS prepare -c ${CONF_OUT} -t ${META_OUT}
#   echo "gemBS will run the following commands:"
#   gemBS --dry-run run
#   cd ${home}
# }

# remove_tmps() {
#   cd ${idir}
#   cd ${MAP_PATH}
#   echo "Removing temp files"
#   rm *.tmp*
#   cd ${home}
# }
