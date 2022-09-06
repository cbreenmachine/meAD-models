#!/bin/bash
all_chr=$(echo "chr"{1..22})
parallel --jobs 4  \
    python impute_missing_compute_PCs.py --chrom {1} ::: $all_chr
#END
