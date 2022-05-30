#!/bin/bash
input=$1
output=$2
while read -r line
do
	fastq-dump $line -O $output
done < $input
