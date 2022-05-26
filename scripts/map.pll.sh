#!/bin/bash
# Sends gemBS map command to several servers, which are hard coded into this command

idir="${1}"

working_dir=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

# echo "calling script from: "${parent_path}
# echo "parallel script housed in "$(pwd)


cd ${idir}
gemBS prepare -c ${CONF_OUT} -t ${META_OUT}


parallel --dry-run \
    -S mastodon-1,mastodon-2,mastodon-3 \
    --joblog ${date}-map.log \
    --nonall --workdir ${working_dir} gemBS map --no-merge

