#!/usr/bin/bash

tmp_file="index.html"
url="https://gg-sg-web.s3-us-west-2.amazonaws.com/"

wget $url -O $tmp_file
files=$(grep -oE "<Key>[^<]*gg_13_5/[^<]*</Key>" $tmp_file | sed 's/<Key>\([^<]*\)<\/Key>/\1/' | grep ".gz$")

for f in $files; do  wget -P ./data/ "$url$f"; done
tar -xvzf ./data/*tar*
gunzip ./data/*.gz
wget -P ./data/ https://zenodo.org/record/5885631/files/sorted_distance_complete.txt
wget -P ./data/ https://zenodo.org/record/5915162/files/mapping_file2.txt
rm $tmp_file

