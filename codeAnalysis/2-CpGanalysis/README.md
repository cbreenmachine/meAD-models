# README.md

Updated extensively December 20, 2022. This is when I split this project into four repos: meAD-fq2bed, meAD-models, meAD-experiments, meAD-paper.

Briefly,
- `0-prepare-DSS-inputs.R` creates the BSseq object that DSS needs, filters based on median coverage, and matches the samplesheet.  
- `1-fit-DSS-models.R` wraps DSS, and exports the fitted $\beta$ values, $p$-values, etc.  
- `2-export-pi-and-p.R`  