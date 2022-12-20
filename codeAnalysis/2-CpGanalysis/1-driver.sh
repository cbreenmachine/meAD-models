#!/bin/bash
# RUNALL.sh
# November 2, 2022
# Does
# New needs: flexibility in DSS modeling script (i.e. fitting models with difference covariates tested)
# Changes: data munging / prep is siloed into 0-prepareDataForDSS.R

njobs=6
all_chroms=$(echo chr{1..22})

# Input directories for data prep step and model fit step rspectively
prep_idirs="../../dataDerived/methylBedImputed-controlLOAD/ ../../dataDerived/methylBedImputed-controlMCILOAD/"
fit_idirs="../../dataDerived/analysis-controlLOAD/ ../../dataDerived/analysis-controlMCILOAD/"
test_vars="diagnostic_group_coded ravlt_long"

date=$(date '+%Y-%m-%d')

# STEP 1: PREPARE DATA (Checks if it's already done in R script)
# parallel --jobs $njobs \
#         --memsuspend 30G \
#         --joblog "logs/${date}.prep.log" \
#         Rscript 0-prepareDataForDSS.R --idir {1} --chr {2} \
#         ::: $prep_idirs ::: $all_chroms

# STEP 2: FIT MODELS (ONLY NEED TO DO ONCE)
parallel --jobs $njobs \
        --joblog "logs/${date}.fit.log" \
        Rscript 1-fitModels.R --idir {1} --chr {2} --test_covariate {3} \
        ::: $fit_idirs ::: $all_chroms ::: $test_vars