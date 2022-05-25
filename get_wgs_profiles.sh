#!/bin/bash
input_dir=$1
out_dir="$input_dir"/profiles
#mode=$3 #precision, recall or parenthesis
[[ -d "$input_dir"/profiles ]] | mkdir $out_dir
for f in `find $input_dir -name "*fastq"` 
	do
		echo $f
		file_name="${f%%.*}"
		bn="${file_name##*/}"
		outfile="$out_dir"/$bn.profile
		echo $outfile
		motus profile -s $f -o $outfile -t 25 -C precision -l 45 -g 1	
done
