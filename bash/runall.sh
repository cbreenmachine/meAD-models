#!/bin/bash

# A nice idiom to break jobs across clusters
for p in $(echo 0{1..9}); do ./extract.pll.sh "../data/pool$p/" 5; done 
for p in $(echo 1{0..9}); do ./extract.pll.sh "../data/pool$p/" 10; done 