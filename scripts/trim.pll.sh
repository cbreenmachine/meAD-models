#!/bin/bash
# Parallelized trimming

#https://stackoverflow.com/questions/24112727/relative-paths-based-on-file-location-instead-of-current-working-directory
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"

idir="${1}"
odir="${2}"

echo ${idir}
    
# Temporary files
> LEFT && > RIGHT

# The following loop is a hacky way to decide if a sample needs trimming
# If there are 4 files starting with 123, then it is assumed that sample
# "completed" trimming. This may not be the case due to incomplete writing,
# storage issues, etc. But as a heuristic works pretty well
# It then curates the files into LEFT (for left / forward reads) 
# and RIGHT (right/reverse) to feed to trim_galore.
# We've stuck a GNU parallel wrapper to speed up

for ff in $(find ${idir} -name *R1.fastq.gz); do 
    un=$(basename ${ff} | cut -c1-3) # grab the 123 from 00-fastq/123_stuff.fastq.gz

    # a successful trim_galore run produces 4 files starting with 123
    # if less than 4 files, needs trimming
    if [[ $(find "${odir}" -name "${un}*" | wc -l) -lt 4 ]] 
    then
        echo "${ff}" >> LEFT; 
        # | sed -e s/val_1/val_2/ 
        echo ${ff} | sed -e s/R1/R2/ >> RIGHT; 
    else 
        echo "${ff}" >> TRIMMED
    fi
    
done

cat TRIMMED|uniq > TRIMMED

head LEFT > LEFT
head RIGHT > RIGHT

#TODO: working dir taken from parent
#TODO: 
# --link creates a mapping between the lines in LEFT and lines in RIGHT 
# (one-to-one like python's zip instead of pairwise combinations)
parallel --link \
    -S nebula-5 \
    --jobs 3 --workdir . \
    --joblog trim.log \
        trim_galore --phred33 --overlap 3 \
        --output_dir ${odir} \
        --dont_gzip --paired {1} {2} :::: LEFT :::: RIGHT