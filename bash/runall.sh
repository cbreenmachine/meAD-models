#!/bin/bash

# A nice idiom to break jobs across clusters

# (gem) environment
for p in $(echo 1{0..9}); do ./map.pll.sh "../data/pool$p/" 3; done 


for p in $(echo {01..20}); do ./extract.pll.sh "../data/pool$p/" 5; done 

# Any environment will do, so long as parallel is installed
./split-by-chrom.pll.sh ../dataDerived/methylBedBySample/ ../dataDerived/methylBedByChrom/ 5
