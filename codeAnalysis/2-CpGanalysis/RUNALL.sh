#!/bin/bash
# RUNALL.sh
# November 2022; Updated December 21.
# New needs: flexibility in DSS modeling script (i.e. fitting models with difference covariates tested)
# Changes: data munging / prep is siloed into 0-prepareDataForDSS.R

njobs=6
all_chroms=$(echo chr{1..22})

# Input directories for data prep step and model fit step rspectively
prep_idirs="../../dataDerived/methylBedImputed-controlLOAD/ ../../dataDerived/methylBedImputed-controlMCILOAD/"
fit_idirs="../../dataDerived/analysis-controlLOAD/ ../../dataDerived/analysis-controlMCILOAD/"
test_vars="diagnostic_group_coded"

date=$(date '+%Y-%m-%d')

# STEP 1: PREPARE DATA 
parallel --jobs $njobs \
        --joblog "logs/${date}.prep.log" \
        Rscript 0-prepare-DSS-inputs.R \
        --idir {1} --chr {2} \
        ::: $prep_idirs ::: $all_chroms

# STEP 2: FIT MODELS (ONLY NEED TO DO ONCE)
parallel --jobs $njobs \
        --joblog "logs/${date}.fit.log" \
        Rscript 1-fit-DSS-models.R \
        --idir {1} --chr {2} --test_covariate {3} \
        ::: $fit_idirs ::: $all_chroms ::: $test_vars