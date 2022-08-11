#!/bin/bash
# Sketch of processing pipeline; usually these steps are run separately

source configuration.txt

./trim.pll.sh ${FASTQ_PATH} ${FASTQ_TRIMMED_PATH}

./extract.pll.sh ${CALLS_PATH} ${EXTRACT_PATH}

./prepare.sh ${idir} ${CONF_OUT} ${META_OUT}

./map.sh ${idir}

find ../data/04-extract/ -type f -name *bed \
    | parallel './split-by-chrom.sh {} ../data/06-counts-by-chrom/'