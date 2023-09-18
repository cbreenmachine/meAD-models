#!/bin/bash
# Handles parallelization of data prep and model fit (DSS)
# See run.config.sh for directory structure
# Last update: September 18, 2023

source run.config.sh; source run.params.sh

# STEP 0: PREPARE DATA 
parallel --jobs $njobs \
        Rscript 0-prepare-DSS-inputs.R \
        --idir ${idir} \
        --odir ${inputs_dir} \
        --chr {1} ::: $all_chroms

# STEP 1: FIT MODELS 
parallel --jobs $njobs \
        Rscript 1-fit-DSS-models.R \
        --idir ${inputs_dir} \
        --odir ${outputs_dir} \
        --test_covariate ${test_var} \
        --chr {} ::: ${all_chroms} 

# END