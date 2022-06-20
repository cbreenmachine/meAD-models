#!/bin/bash

root_dir="$1"
conf_file="$2"
meta_file="$3"

cd ${root_dir}

echo "barcode,dataset,end1,end2" > ${meta_file}
# 216,216,216_val_1.fq,216_val_2.fq
find 01-fastq-trimmed/ -type f -name *val_1.fq -printf "%f\n" \
    | cut -d _ -f1 \
    | sed 's/\([0-9][0-9][0-9]\)/\1,\1,\1_R1_val_1\.fq,\1_R2_val_2.fq/' >> ${meta_file}

gemBS prepare -c ${conf_file} -t ${meta_file}