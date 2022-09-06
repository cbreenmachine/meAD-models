# root/bash/

This directory mostly contains bash scripts used to manipulate sequencing files: fastq, sam/bam, vcf/bcf. Scripts are meant to run parallel instances if requested by analyst. These parallelized scripts end in `.pll.sh` for parallel and use GNU parallel. It's recommended to install GNU parallel and add to you path outside of conda package management. You'll also need a working copy of samtools and bcftools. Conda often messes up the shared libraries between these softwares (see for example (this github issue)[https://github.com/bioconda/bioconda-recipes/issues/12100#issuecomment-461466897]).

## Mapping
`./map.pll.sh ../data/pool01/ 3` runs three parallel mapping instances of gem-mapper. Each instance needs 30GB of memory and 4 CPUs (can be changed in the wrapper script). On Mastodon-x, can run 10 or so parallel instances. Produces files like `../data/pool01/100.unsorted.sam`.

## Soring
`./sort.pll.sh ../data/pool01/ 3` sorts and binarizes sam files from mapping. We don't pipe directly to sort because gem-mapper fails quite regularly, and piping to `samtools sort` often delays the error raise. Produces files like `../data/pool01/100.bam`.

## Calling

Uses bs_call and runs like `./call.pll.sh ../data/pool01/ 3`. Produces files like `../data/pool01/100.bcf`. May parallelize by chromosome so you don't need to wait a long time for one file.

## Extaction 

Run in (wgbs). This is the final step and takes the most pertinent information (i.e. methylated versus not-methylated) and outputs a bed file.

To prevent search-path issues, `extract_from_bcf.py` is stored in this (bash) directory. This is mostly for clarity, but may change as the project evolves. As many will recognize immediately, increasing the chunksize like `df.to_csv(args.ofile, sep='\t', chunksize=30000000)` is extremely important to speed. The `bcftools query` command takes about one hour to write to a text file (assumption is that writing is making about something like half of this). The python script with the chunksize modification runs in <8 minutes (as opposed to 5 hours on default).

For some mystical reason, the following code block (note the exclusion filter) works.
```
bcftools query \
    --regions chr1 \
    --format '%CHROM\t%POS\t[%CS]\t%REF\t[%MC8{4}]\t[%MC8{5}]\t[%MC8{6}]\t[%MC8{7}]\n' \
    --exclude 'CG="N" | CG="H" | CG="?"' ../data/pool20/499.bcf 
```
while its equivalent (*inclusion* filter) does not.

```
bcftools query \
    --regions chr1 \
    --format '%CHROM\t%POS\t[%CS]\t%REF\t[%MC8{4}]\t[%MC8{5}]\t[%MC8{6}]\t[%MC8{7}]\n' \
    --exclude 'DP<=2 | DP>=50 | CG="N" | CG="H" | CG="?"' "${ifile}" 
```

It takes 20-25 mins to query (and print to nothing), and about 8 mins to extract.

## SNP-extraction

TBD, but we'll need to pull out SNP information. 

## Driver 

A driver script `RUNALL.sh` shows how one *could* run these scripts in order. In reality, there was a lot of "mapping and calling on this subset, sorting on this subset" type processing.

# Note on 0- and 1-based coordinate systems

| Step | Name | Input Type | Ouptut Type | Coordinate of Output |
|------|------|------------|-------------|----------------------|
| 1 | Mapping | fastq | sam |