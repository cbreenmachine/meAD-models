#!/bin/bash
# November 7, 2022
# Does: runs pi calculation (adjusted response from DSS) on 
# all control-LOAD CpGs, one chromosome at a time.
# The internals of 2-computePis.R are parallelized, so no need to run GNU parallel

all_chroms=$(echo chr{1..22})

# Input directories for data prep step and model fit step rspectively
idir="../../dataDerived/analysis-controlLOAD/test-diagnostic-group-coded/"

for c in $all_chroms;
do 
    Rscript 2-computePis.R --idir ${idir} --chr ${c}
done