#!/bin/bash

# Plot PCs for just array
# Plot PCs for all

Rscript plot-PCs.R --idir ../../data/prin-comps/ --odir ./figs/
Rscript plot-PCs.R --idir ../../data/prin-comps-array/ --odir ./figs-array/