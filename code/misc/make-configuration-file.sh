
ifile=../data/master-samplesheet.csv

ofile=../data/meta.csv

echo "barcode,dataset,end1,end2" > ${ofile}

for s in $(tail -n +2 "${ifile}" | awk -F',' '$2 != "NA" {print $2}');
do
echo "${s},${s},${s}_val_1.fq,${s}_val_2.fq" >> ${ofile}

done