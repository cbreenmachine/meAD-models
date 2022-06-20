#!/bin/bash
all_chr=$(echo "chr"{1..22})
parallel --jobs 4  \
    python impute_missing_compute_PCs.py \
        --idir ../data/06-counts-by-chrom/ \
        --samples_file ../data/07-counts-by-chrom-imputed-subset1/master-42-subsampled-1.csv \
        --odir ../data/07-counts-by-chrom-imputed-subset1/ \
        --chrom {1} ::: $all_chr
#END