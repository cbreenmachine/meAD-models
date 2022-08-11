#!/bin/bash

declare -A files_to_download

#files_to_download[Alisch_Pool16_group1.202216.tar.bz2]="https://posting.biotech.illinois.edu/posting/alisch/Alisch_Pool16_group1.202216.tar.bz2?x-email=alisch%40neurosurgery.wisc.edu&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=posting%2F20220106%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20220106T183701Z&X-Amz-Expires=604800&X-Amz-SignedHeaders=host&X-Amz-Signature=0b221a740a44e7d3979bbd4005f65a1b4300fcc7b2e87839328829df1ec44c51"
#files_to_download[Alisch_Pool16_group2.202216.tar.bz2]="https://posting.biotech.illinois.edu/posting/alisch/Alisch_Pool16_group2.202216.tar.bz2?x-email=alisch%40neurosurgery.wisc.edu&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=posting%2F20220106%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20220106T183701Z&X-Amz-Expires=604800&X-Amz-SignedHeaders=host&X-Amz-Signature=48fa54ef1cf72c606b298a335dd0d90c23f12bbded84e1e013841e50a2f49825"
#files_to_download[Alisch_Pool18_group1.202215.tar.bz2]="https://posting.biotech.illinois.edu/posting/alisch/Alisch_Pool18_group1.202215.tar.bz2?x-email=alisch%40neurosurgery.wisc.edu&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=posting%2F20220106%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20220106T183642Z&X-Amz-Expires=604800&X-Amz-SignedHeaders=host&X-Amz-Signature=d9874ff84de73793349b37cc0e2a09d154c321664fe4727663a78db76fb2d2c9"
#files_to_download[Alisch_Pool18_group2.202215.tar.bz2]="https://posting.biotech.illinois.edu/posting/alisch/Alisch_Pool18_group2.202215.tar.bz2?x-email=alisch%40neurosurgery.wisc.edu&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=posting%2F20220106%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20220106T183642Z&X-Amz-Expires=604800&X-Amz-SignedHeaders=host&X-Amz-Signature=59ed82435b8affb6137230d9059c69f45aedd1d32935c82a5af6f3adf6708442"
files_to_download[Alisch_Pool17_group1.2022111.tar.bz2]="https://posting.biotech.illinois.edu/posting/alisch/Alisch_Pool17_group1.2022111.tar.bz2?x-email=alisch%40neurosurgery.wisc.edu&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=posting%2F20220114%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20220114T025301Z&X-Amz-Expires=604800&X-Amz-SignedHeaders=host&X-Amz-Signature=0ca2246490d3a587fec9e5ea902c76bb4d08113d57916cce2e12984c6f5ebe2d"
#files_to_download[Alisch_Pool17_group2.2022111.tar.bz2]="https://posting.biotech.illinois.edu/posting/alisch/Alisch_Pool17_group2.2022111.tar.bz2?x-email=alisch%40neurosurgery.wisc.edu&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=posting%2F20220114%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20220114T023923Z&X-Amz-Expires=604800&X-Amz-SignedHeaders=host&X-Amz-Signature=bbdfe1c1f4b7d39091c595a0a06ce8b6a741d133e17e8e96831877aaaac7e6c5"
#files_to_download[Alisch_Pool19.2022114.tar.bz2]="https://posting.biotech.illinois.edu/posting/alisch/Alisch_Pool19.2022114.tar.bz2?x-email=alisch%40neurosurgery.wisc.edu&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=posting%2F20220115%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20220115T181011Z&X-Amz-Expires=604800&X-Amz-SignedHeaders=host&X-Amz-Signature=c1da4283f7ec78598858ca7de1a148da92f2b9a7e6d96b4db9ba90368c792c56"
#files_to_download[Alisch_Pool20.2022114.tar.bz2]="https://posting.biotech.illinois.edu/posting/alisch/Alisch_Pool20.2022114.tar.bz2?x-email=alisch%40neurosurgery.wisc.edu&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=posting%2F20220115%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20220115T172644Z&X-Amz-Expires=604800&X-Amz-SignedHeaders=host&X-Amz-Signature=775e88e1487a21254318f000af4c106f998d97bce4ecd0195eba40dc9237932a"

function download_retry() {
    while true; do
        echo ${1}
        wget -T 15 --continue -O "${1}" "${2}" && break
    done
}

for c in "${!files_to_download[@]}"; do
    echo ${c}
    download_retry "${c}" "${files_to_download[$c]}"
done
