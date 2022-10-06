#!/bin/bash
# Paths handled in python directory
# each instance needs ~10GB / 20GB virtual
all_chr=$(echo "chr"{1..22})
njobs=11

# No MCI
parallel --jobs $njobs --joblog logs/impute.controlLOAD.log \
    python impute_missing_compute_PCs.py --chrom {1} ::: $all_chr

# With MCI
parallel --jobs $njobs --joblog logs/impute.controlMCILOAD.log \
    python impute_missing_compute_PCs.py --include_MCI \
    --chrom {1} ::: $all_chr
#END