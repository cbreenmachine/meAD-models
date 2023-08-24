#!/bin/bash

njobs=6
all_chroms=$(echo chr{1..22})

# Should be 
# Inputs/
# test-diagnostic-group-coded/
#       Outputs/
#       Pis/
#       Summaries/

idir="../DataRaw/MCovRawSplit/"
root_odir="../DataDerived/ControlLOAD/"

inputs_dir="${root_odir}/Inputs/"
test_dir="${root_odir}/test-diagnostic-group-coded/"

mkdir -p $test_dir

# Nested under test-covariate-of-interest
outputs_dir="${test_dir}/Outputs/"
pis_dir="${test_dir}/Pis/"
summaries_dir="${test_dir}/Summaries/"

mkdir -p $outputs_dir
mkdir -p $pis_dir
mkdir -p $summaries_dir

# This needs to be in the format found in data frame
test_var="diagnostic_group_coded"
date=$(date '+%Y-%m-%d')
# END