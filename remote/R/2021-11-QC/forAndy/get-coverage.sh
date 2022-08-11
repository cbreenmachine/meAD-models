ifile="../../../data/batch01/pool01-group02/03-calls/110.bcf"

bcftools query \
    -f '[ %DP]\n' \
    --include "DP >= 2 & DP <= 80 && CG = 'Y'" \
    ${ifile} > 110.cpg.cov