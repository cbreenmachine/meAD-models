#!/bin/bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.chr_patch_hapl_scaff.basic.annotation.gtf.gz
gunzip gencode.v40.chr_patch_hapl_scaff.basic.annotation.gtf.gz

mv gencode.v40.chr_patch_hapl_scaff.basic.annotation.gtf data/