#!/bin/bash
#https://broadinstitute.github.io/picard/explain-flags.html
#2: mapped as pair
samtools view -h -f 2 261.trimmed.sorted.bam \
    | bs_call --reference ../../../reference/GENCODE/h38_no_alt.fa -L4 > 262.trimmed.bcf