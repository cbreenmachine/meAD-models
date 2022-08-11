# https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz
# https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.chr_patch_hapl_scaff.basic.annotation.gtf.gz
# awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id "";"; }' \
#     gencode.v40.chr_patch_hapl_scaff.basic.annotation.gtf \
#     | gtf2bed > gencode.hg38.bed

bedparse gtf2bed gencode.v40.chr_patch_hapl_scaff.basic.annotation.gtf > gencode.hg38.bed