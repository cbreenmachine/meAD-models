#!/bin/bash
# Sends gemBS map command to several servers, which are hard coded into this command

date=$(date '+%Y-%m-%d')
cd ${1}

parallel -S mastodon-1,mastodon-2 \
    --joblog ${date}-map.log \
    --nonall --workdir . \
    gemBS --dry-run map --no-merge
