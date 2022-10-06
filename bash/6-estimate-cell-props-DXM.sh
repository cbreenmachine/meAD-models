#!/bin/bash
# ./6-estimate-cell-props-DXM.sh ../dataDerived/methylBedBySample/ ../dataDerived/methylBedForDXM/ 10
idir="${1}"
odir="${2}"
njobs="${3}"

# Where we stroe all the estimates once they're computed
ofile="${odir}/bloodCellPropsDXM.csv"

find "${idir}" -name *bed |sort > dxm_inputs
find "${idir}" -name "*bed" -printf '%f\n' | sed s/bed/dxm\.bed/ | sed -e "s@^@$odir@g" > dxm_outputs
find "${idir}" -name "*bed" -printf '%f\n' | sed s/\.bed// | sed -e "s@^@$odir@g" > dxm_prev_outputs

# Phase 1: get data into format that works with DXM. It needs a BED with specific columns
format_wrapper() {
    ifile="$1"
    ofile="$2"

    if [ -f "$ofile" ]
    then
        echo "$ofile exists; delete if you want to run formating"
    else
        python format_bed_for_dxm.py --ifile "$1" --ofile "$2"
    fi
}

export -f format_wrapper
parallel --jobs "${njobs}" --link \
    format_wrapper {1} {2} :::: dxm_inputs :::: dxm_outputs
    

estimate_wrapper() {
    ifile="$1"
    ofile="$2"

    if [ -f "$ofile" ]
    then
        echo "$ofile exists; delete if you want to estimate cell composition"
    else
        dxm_estimateFracs -i $ifile -k 6 -o $ofile
    fi
}

# Phase 2: estimate proportions
parallel --jobs "${njobs}" --link \
    dxm_estimateFracs -i {1} -k 6 -o {2} :::: dxm_outputs :::: dxm_prev_outputs

rm dxm_inputs dxm_outputs dxm_prev_outputs

# Phase 3: Collate into one file
echo "sample_id,type1,type2,type3,type4,type5,type6" > $ofile

for ff in $(ls ${odir}*txt)
do
    sample_id=$(basename $ff | cut -d "_" -f1)
    nums=$(cat $ff | tr '\n' ',' | rev | cut -c2- | rev) 

    echo "$sample_id,$nums" >> $ofile
done
#END