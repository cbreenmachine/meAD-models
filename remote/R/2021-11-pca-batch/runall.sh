#!/bin/bash

for ii in $(echo {1..22});
do
    Rscript plot-PCs.R --chrom "chrom${ii}"
done
