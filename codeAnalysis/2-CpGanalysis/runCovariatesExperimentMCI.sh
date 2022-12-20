#! /bin/bash

# PCs vs no PCs
# cell type composition vs no cell type composition
# 

njobs=11


#CONSTANTS
all_chroms=$(echo chr{1..22})

idir="../../dataDerived/methylBedImputed-controlMCILOAD/"
models_root="../../dataDerived/models"
dmps_root="../../dataSummaries/DMPs-"
config="controlMCILOAD"

baseline_covariates="diagnostic_group_coded,bmi,sex,age_at_visit"


parallel --jobs $njobs \
    --joblog "logs/model.${config}.pc0-ct0.log" \
    Rscript 1-fitModels.R \
        --idir "${idir}" \
        --odir "${models_root}${config}-pc0-ct0/" \
        --config "${config}" \
        --covariates "${baseline_covariates}" \
        --num_pcs 0 \
        --chr {1} ::: $all_chroms


parallel --jobs $njobs \
    --joblog "logs/model.${config}.pc2-ct0.log" \
    Rscript 1-fitModels.R \
        --idir "${idir}" \
        --odir "${models_root}${config}-pc2-ct0/" \
        --config "${config}" \
        --covariates "${baseline_covariates}" \
        --num_pcs 2 \
        --chr {1} ::: $all_chroms


parallel --jobs $njobs \
    --joblog "logs/model.${config}.pc0-ct5.log" \
    Rscript 1-fitModels.R \
        --idir "${idir}" \
        --odir "${models_root}${config}-pc0-ct5/" \
        --config "${config}" \
        --covariates "${baseline_covariates},type1,type2,type3,type4,type5" \
        --num_pcs 0 \
        --chr {1} ::: $all_chroms


parallel --jobs $njobs \
    --joblog "logs/model.${config}.pc2-ct5.log" \
    Rscript 1-fitModels.R \
        --idir "${idir}" \
        --odir "${models_root}${config}-pc2-ct5/" \
        --config "${config}" \
        --covariates "${baseline_covariates},type1,type2,type3,type4,type5" \
        --num_pcs 2 \
        --chr {1} ::: $all_chroms



# Export DMPs, make figs
idir=../../dataDerived/models-${config}
odir=../../dataSummaries/${config}-
fig_dir=../../figs/pvals-${config}-
s_range="pc0-ct0 pc2-ct0 pc0-ct5 pc2-ct5"

for s in $(echo $s_range); do
    desc=$(echo $s | sed 's/pc\([0-2]\)-ct\([0-5]\)/\1 PCs and \2 cell types/')
    Rscript 2-exportDMPsAsBED.R \
        --idir "${idir}-${s}/diagnostic_group_coded/" \
        --odir "${odir}${s}" \
        --fig_dir "${fig_dir}${s}" \
        --title "${desc}" 
done