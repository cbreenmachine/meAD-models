
ofile=bloodCellPropsDXM.csv

echo "sample_id,type1,type2,type3,type4,type5,type6" > $ofile

for ff in ../dataDerived/methylBedForDXM/*
do
    sample_id=$(basename $ff | cut -d "_" -f1)
    nums=$(cat $ff | tr '\n' ',' | rev | cut -c2- | rev) 

    echo "$sample_id,$nums" >> $ofile
done
