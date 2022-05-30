#!/bin/bash
input_dir=$1
manifest="$input_dir/16s_reads/manifest.txt"
echo "sample-id,absolute-filepath,direction" > $manifest
for f in `find $input_dir/16s_reads -name "*fastq"`
do
	file_name="${f%%.*}"
        id="${file_name##*/}"
	printf '%s,%s,%s\n' "$id" "$PWD/$f" "forward" >> $manifest
done
