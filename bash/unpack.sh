ls tarballs/*.bz2 | cut -d "_" -f2 | sed -e s/P/p/ | sed 's/$/\/00-fastq\//' > DIR
ls tarballs/*.bz2 > FILE

parallel --jobs 5 --link tar xvkf {1} --driectory {2} :::: FILE :::: DIR