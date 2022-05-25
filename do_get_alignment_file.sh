#!/bin/bash
###convert all read files in a particular experiment into .sam files 
###useage: bash do_get_alignment_files.sh experiment_dir
###requirements: bowtie2 installed and working, database present (here, ./bowtie2/genome2)

input=$1 #experiment directory
reads_dir="$input/reads"
out_dir="$input/bt2out"
db="./genomes/wol_mini"

#create out_dir if doesn't exist
[ -d $out_dir ] || mkdir $out_dir

for rf in $reads_dir/*fastq
do
	file_name="${rf%%.*}"
        bn="${file_name##*/}"
	echo $bn
	bowtie2 --very-sensitive -p 8 -x $db -q $rf -S "$out_dir/$bn.sam"
done
