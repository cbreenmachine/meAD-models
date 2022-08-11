#!/bin/bash
# Sketch of processing pipeline; usually these steps are run separately

./process trim_all

./process extract_all
./process split_by_chrom