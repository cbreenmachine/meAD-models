#!/bin/bash

effect_range="0 0.1 0.25 0.5 0.75"
upstream_range="0"

parallel --jobs 5 \
    Rscript 3-callDMGenes.R \
    --effect_cut {1} --upstream {2} \
    ::: $effect_range ::: $upstream_range
