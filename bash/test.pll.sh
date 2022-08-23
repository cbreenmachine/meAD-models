#!/bin/bash

# my_func(){
#     conda activate gem3
#     conda info --envs
#     echo $HOSTNAME
# }
# issue is the environmetn
# export -f my_func
~/bin/parallel -S nebula-1,nebula-2 \
    --workdir /z/Comp/kelesgroup/loadmethyseq/methylome-diff-analysis/bash/ \
    ./my_func.sh ::: 1 2 3 4