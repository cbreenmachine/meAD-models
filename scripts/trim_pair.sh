#!/bin/bash

left="${1}"
right="${2}"

odir="${3}"

trim_galore --phred33 --overlap 3 \
        --output_dir ${odir} \
        --dont_gzip --paired ${left} ${right}