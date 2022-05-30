#!/bin/bash
input_dir=$1
for dir in $input_dir/*
do
	abund_files=$dir/wgs_abundance_files
	echo $abund_files
	out_dir=$dir/wgs_reads
	bash get_wgs_libraries.sh $abund_files $out_dir &
done
