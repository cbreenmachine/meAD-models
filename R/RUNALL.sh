#!/bin/bash
# RUNALL.sh
# November 2022; Updated December 21.
# New needs: flexibility in DSS modeling script (i.e. fitting models with difference covariates tested)
# Changes: data munging / prep is siloed into 0-prepareDataForDSS.R

njobs=4
all_chroms=$(echo chr{1..22})

dir1="../DataRaw/MCovRawSplit/"
dir2="../DataDerived/ControlLOAD/Inputs/"
test_var_range="diagnostic_group_coded"

date=$(date '+%Y-%m-%d')

# STEP 0: PREPARE DATA 
parallel --jobs $njobs \
        --joblog "logs/${date}.prep.log" \
        Rscript 0-prepare-DSS-inputs.R \
        --idir {1} --odir {2} --chr {3} \
        ::: $dir1 ::: $dir2 ::: $all_chroms

# STEP 1: FIT MODELS 
parallel --jobs $njobs \
        --joblog "logs/${date}.fit.log" \
        Rscript 1-fit-DSS-models.R \
        --idir {1} --chr {2} --test_covariate {3} \
        ::: ${dir2} ::: ${all_chroms} ::: ${test_var_range}

# STEP 2: Extract p-values and pi (corrected methylation percentages)
wrapper() {
        path="${1}"
        v="${2}"

        # String maniupulation to match output directories
        last_dir=$(echo ${v} | sed 's/\_/\-/g')
        my_idir="${path}/test-${last_dir}"

        # Computationally expensive bit
        Rscript 2-extract-pi-and-pvals.R --idir ${my_idir} 
}

export -f wrapper
parallel --jobs 1 \
        --joblog "logs/${date}.extract.log" \
        wrapper ${1} ${2} \
        ::: ${dir2} ::: ${test_var_range} 



p1=${dir2}"/test-"(sed s/${test_var}_/-/)


# END