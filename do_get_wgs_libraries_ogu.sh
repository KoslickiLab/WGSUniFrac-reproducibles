#!/bin/bash
dir=$1 #the e2r... directory
db=data/concat.fna
#generates wgs reads for each experiment
#[ -d $dir/reads ] | mkdir $dir/reads

N=10
for f in $dir/abundance_files/*
do
	file_name="${f%%.*}"
	bn="${file_name##*/}"
	#echo $bn
	#profile=$bn-reads.profile
	reads=$bn-reads.fastq
	output_dir=$dir/reads
	#echo $profile
	#if [ -f $dir/profile_precision/$profile ] #if profile exists
	#then
	#	echo "profile exists"
	#	rm $dir/wgs_reads/$reads || true #remove reads 
	#elif  [ -f $dir/wgs_reads/$reads ] && [ "$(cat "$dir"/wgs_reads/"$reads" | wc -l)" -eq 4000000 ] 
	#then
	#	echo "reads ready"
	#echo "generate reads"
	grinder -rf $db -tr 1000000 -af $f -md poly4 3e-3 3.3e-8 -mr 80 20 -ql 50 10 -bn $bn -od $output_dir -fq 1 -rd 150 &
	(( ++count % N == 0)) && wait
	
done
wait
echo "all done"
