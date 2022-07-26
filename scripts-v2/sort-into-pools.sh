#!/bin/bash


for pool in $(echo pool0{1..9}); do
    for num in $(grep ${pool} processed.txt | cut -c1-3); do
        find ../data/03-calls/ -name "${num}*" -type f -exec mv {} ../data/${pool} \;
    done
done
