# ifile=../data/sample-info/master-samplesheet.csv
idir=../data/02-mapping
ofile=../data/meta.csv
trimmed_path=../data/01-fastq-trimmed/

echo "barcode,dataset,end1,end2" > ${ofile}

# for s in $(tail -n +2 "${ifile}" | awk -F',' '$1 != "NA" {print $1}');
# for s in $(find ${map_dir} -name *bam -printf "%f\n" | sed s/\.bam//)
for s in $(find ${trimmed_path} -name *.fq -printf "%f\n" |cut -d _ -f1 |sort |uniq)
do
    echo "${s},${s},${s}_R1_val_1.fq,${s}_R2_val_2.fq" >> ${ofile}
done


# for s in $(tail -n +2 "${ofile}");
# do

# left_file="${trimmed_path}"/$(echo "${s}" | cut -d , -f3)
# right_file="${trimmed_path}"/$(echo "${s}" | cut -d , -f4)

#     if [ ! -f "${left_file}" ]
#     then
#         > "${left_file}"
#         > "${right_file}"
#     fi
# done
