#! /bin/bash

# PCs vs no PCs
# cell type composition vs no cell type composition
# 

njobs=11


#CONSTANTS
all_chroms=$(echo chr{1..22})

idir="../../dataDerived/methylBedImputed-controlLOAD/"
models_root="../../dataDerived/models"
dmps_root="../../dataSummaries/DMPs-"
config="controlLOAD"

baseline_covariates="diagnostic_group_coded,bmi,sex,age_at_visit"


parallel --jobs $njobs \
    --joblog "logs/model.pc0-ct0.log" \
    Rscript 1-fitModels.R \
        --idir "${idir}" \
        --odir "${models_root}-pc0-ct0/" \
        --config "controlLOAD" \
        --covariates "${baseline_covariates}" \
        --num_pcs 0 \
        --chr {1} ::: $all_chroms


parallel --jobs $njobs \
    --joblog "logs/model.pc2-ct0.log" \
    Rscript 1-fitModels.R \
        --idir "${idir}" \
        --odir "${models_root}-pc2-ct0/" \
        --config "controlLOAD" \
        --covariates "${baseline_covariates}" \
        --num_pcs 2 \
        --chr {1} ::: $all_chroms


parallel --jobs $njobs \
    --joblog "logs/model.pc0-ct5.log" \
    Rscript 1-fitModels.R \
        --idir "${idir}" \
        --odir "${models_root}-pc0-ct5/" \
        --config "controlLOAD" \
        --covariates "${baseline_covariates},type1,type2,type3,type4,type5" \
        --num_pcs 0 \
        --chr {1} ::: $all_chroms


parallel --jobs $njobs \
    --joblog "logs/model.pc2-ct5.log" \
    Rscript 1-fitModels.R \
        --idir "${idir}" \
        --odir "${models_root}-pc2-ct5/" \
        --config "controlLOAD" \
        --covariates "${baseline_covariates},type1,type2,type3,type4,type5" \
        --num_pcs 2 \
        --chr {1} ::: $all_chroms



# Export DMPs, make figs
idir=../../dataDerived/models-
odir=../../dataSummaries/
fig_dir=../../figs/
s_range="pc0-ct0 pc2-ct0 pc0-ct5 pc2-ct5"

for s in $(echo $s_range); do
    desc=$(echo $s | sed 's/pc\([0-2]\)-ct\([0-5]\)/\1 PCs and \2 cell types/')
    Rscript 2-exportDMPsAsBED.R \
        --idir "${idir}${config}-${s}/diagnostic_group_coded/" \
        --odir "${odir}dmps-${config}-${s}" \
        --fig_dir "${fig_dir}pvals-${config}-${s}" \
        --title "${desc}" 
done