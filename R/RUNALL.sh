#!/bin/bash
# RUNALL.sh
# November 2022; Updated December 21.
# New needs: flexibility in DSS modeling script (i.e. fitting models with difference covariates tested)
# Changes: data munging / prep is siloed into 0-prepareDataForDSS.R

njobs=4
all_chroms=$(echo chr{1..22})

idir_range="../DataRaw/methylBedImputed-controlLOAD/"
odir_range="../DataDerived/controlLOAD/"
test_var_range="diagnostic_group_coded"

date=$(date '+%Y-%m-%d')

# STEP 0: PREPARE DATA 
# parallel --jobs $njobs \
#         --joblog "logs/${date}.prep.log" \
#         Rscript 0-prepare-DSS-inputs.R \
#         --idir {1} --chr {2} \
#         ::: $prep_idirs ::: $all_chroms

# STEP 1: FIT MODELS 
# parallel --jobs $njobs \
#         --joblog "logs/${date}.fit.log" \
#         Rscript 1-fit-DSS-models.R \
#         --idir {1} --chr {2} --test_covariate {3} \
#         ::: $fit_idirs ::: $all_chroms ::: $test_vars

# STEP 2: Extract p-values and pi (corrected methylation percentages)
Rscript 2-extract-pi-and-pvals.R
        --idir {1} 