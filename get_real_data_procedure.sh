#!/bin/bash
input=$1
#log=$2
profile_dir=results/real/profiles
mkdir -p $profile_dir
while read -r line
do
	#echo $line >> $log
	mkdir $PWD/SRA/$line #make directory 
	fastq-dump --split-files $line -O $PWD/SRA/$line #download file
	if [ "$(ls $PWD/SRA/$line | wc -l)" -ne 2 ] #if downloaded does not contain 2 files 
	then 
		echo "not paired. deleted" >> $log
		rm -r $PWD/SRA/$line #delete everything downloaded and continue with the next
	else
		#the 2 files are in the form $line_1.fastq $line_2.fastq
		fastp -i "$PWD/SRA/$line/$line"_1.fastq -I "$PWD/SRA/$line/$line"_2.fastq -o "$PWD/SRA/$line/$line"_1_qc.fastq -O "$PWD/SRA/$line/$line"_2_qc.fastq
		if [[ "$(wc -l "$PWD/SRA/$line/$line"_1_qc.fastq | awk '{print $1}')" -lt 10000 || "$(wc -l "$PWD/SRA/$line/$line"_2_qc.fastq | awk '{print $1}')" -lt 10000 ]] #if after qc too few are left
			then 
				echo "too few after qc. deleted." >> $log
				rm -r $PWD/SRA/$line #delete everything 
		else						
			echo "$(head "$PWD/SRA/$line/$line"_1_qc.fastq)" >> $log 
			motus profile -f "$PWD/SRA/$line/$line"_1_qc.fastq -r "$PWD/SRA/$line/$line"_2_qc.fastq -o "$profile_dir"/"$line".profile -C recall -l 45 -t 20
			#echo "profile generated" >> $log
			rm -r $PWD/SRA/$line #remove after everything is done		
		fi
	fi
done < $input
