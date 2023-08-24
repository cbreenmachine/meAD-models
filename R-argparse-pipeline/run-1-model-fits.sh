#!/bin/bash
# Handles parallelization of data prep and model fit (DSS)
# See run.config.sh for directory structure
# Last update: January 23, 2023

source run.config.sh; source run.params.sh

# STEP 0: PREPARE DATA 
parallel --jobs $njobs \
        --joblog "logs/${date}.prep.log" \
        Rscript 0-prepare-DSS-inputs.R \
        --idir ${idir} \
        --odir ${inputs_dir} \
        --chr {1} ::: $all_chroms

# STEP 1: FIT MODELS 
parallel --jobs $njobs \
        --joblog "logs/${date}.fit.log" \
        Rscript 1-fit-DSS-models.R \
        --idir ${inputs_dir} \
        --odir ${outputs_dir} \
        --baseline_covariates="bmi,sex,age_at_visit,PC1,PC2" \
        --test_covariate ${test_var} \
        --chr {} ::: ${all_chroms} 

# END