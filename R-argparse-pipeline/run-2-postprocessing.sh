#!/bin/bash
# After model fits, all the DMP, DMR calling
# Updated January 24, 2023

source run.config.sh; source run.params.sh

Rscript 2-extract-pi-and-pvals.R \
    --idir ${outputs_dir} \
    --odir_pis ${pis_dir} \
    --odir_summaries ${summaries_dir}

Rscript 3-call-positions.R \
    --ifile "${summaries_dir}/pvals.bed" \
    --ofile "${summaries_dir}/DMPs.bed" \
    --ofile_dmps "${summaries_dir}/DMPs-high-effect.bed"

Rscript 4-call-regions.R \
    --ifile "${summaries_dir}/pvals.bed" \
    --odir "${summaries_dir}"

Rscript 5-subset-pi-by-dmr.R \
    --idir "${pis_dir}" \
    --dmr_file "${summaries_dir}/DMRegions.bed"

Rscript 6-subset-M-Cov-by-dmr.R \
    --idir "${inputs_dir}" \
    --dmr_file "${summaries_dir}/DMRegions.bed"
# END