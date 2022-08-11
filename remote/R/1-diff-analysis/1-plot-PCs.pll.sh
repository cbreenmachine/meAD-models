#!/bin/bash

echo chr{1..22} | tr ' ' '\n' | parallel --jobs 11 Rscript 1-plot-PCs.R --chr {1}
