#!/bin/bash
find ../data/04-extract/ -type f -name *bed \
    | parallel './split-by-chrom.sh {} ../data/06-counts-by-chrom/'