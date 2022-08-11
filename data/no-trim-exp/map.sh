#!/bin/bash
# gem-mapper -I ../../../reference/GENCODE/gembs-index/h38_no_alt.BS.gem \ 
#     --i1 01-fastq/261_R1_val_1.fq --i2 01-fastq/261_R2_val_2.fq --paired-end-alignment \
#     --bisulfite-conversion inferred-C2T-G2A \
#     --report-file 02-mapping/261_trimmed.json \
#     --sam-read-group-header @RG\tID:261_trimmed\tSM:\tBC:261_trimmed\tPU:261_trimmed \
#     | /ua/cebreen/lib/gemBS/bin/read_filter ../../../reference/GENCODE/gembs-index/h38_no_alt.gemBS.contig_md5 \
#     | /ua/cebreen/lib/gemBS/bin/samtools sort -o 02-mapping/261_trimmed.bam -T 02-mapping --write-index -


map_wrapper(){

    # Relative to where GNU parallel running
    ref_path=../../../reference/GENCODE/gembs-index/h38_no_alt.BS.gem

    # Figure out path names here...

    gem-mapper -I "${ref_path}" \
        --report-file 261.json \
        --i1 261_R1.fastq.gz --i2 261_R2.fastq.gz -o 261.unsorted.sam
}