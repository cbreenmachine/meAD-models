#!/bin/bash

Rscript 4-annotateToGenesGivenBED.R \
    --ifile=../../dataSummaries/analysis-controlLOAD/test-dg/dmrs.all.bed \
    --odir=../../figs/2022-paper/dmrGO/ 


# Rscript annotateToGenesGivenBED.R \
#     --ifile=../../dataSummaries/analysis-controlLOAD/test-dg/dmrs.all.bed \
#     --odir=../../dataSummaries/analysis-controlLOAD/test-dg/dmrs.all/ \

# Rscript annotateToGenesGivenBED.R \
#     --ifile=../../dataSummaries/analysis-controlLOAD/test-dg/dmrs.filt.bed \
#     --odir=../../dataSummaries/analysis-controlLOAD/test-dg/dmrs.filt/ \

# Positions
Rscript annotateToGenesGivenBED.R \
    --ifile=../../dataSummaries/analysis-controlLOAD/test-dg/dmps.all.bed \
    --odir=../../dataSummaries/analysis-controlLOAD/test-dg/dmps.all/ \
    --run_pval_correction


Rscript annotateToGenesGivenBED.R \
    --ifile=../../dataSummaries/analysis-controlLOAD/test-dg/dmps.filt.bed \
    --odir=../../dataSummaries/analysis-controlLOAD/test-dg/dmps.filt/ \
    --run_pval_correction