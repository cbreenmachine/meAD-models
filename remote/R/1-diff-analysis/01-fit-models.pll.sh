#!/bin/bash

echo chr{1..22} | tr ' ' '\n' | parallel --jobs 11 Rscript 01-fit-models.R --chr {1}
