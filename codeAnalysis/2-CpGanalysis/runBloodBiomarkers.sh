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
test_vars="a_beta_42_40_bin p_tau_abeta42_bin p_tau_bin"

date=$(date '+%Y-%m-%d')

# STEP 2: FIT MODELS (ONLY NEED TO DO ONCE)
parallel --jobs $njobs \
        --joblog "logs/${date}.fit.log" \
        Rscript 1-fitModels.R --idir {1} --chr {2} --test_covariate {3} \
        ::: $fit_idirs ::: $all_chroms ::: $test_vars