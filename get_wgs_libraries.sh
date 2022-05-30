#!/bin/bash
input_dir=$1 #simulation_batch_....
source /data/wjw5274/anaconda3/etc/profile.d/conda.sh
for dir in $input_dir/* 
do
	bash ./do_get_wgs_libraries.sh $dir </dev/null 
	#conda activate wgsunifrac
	for f in $dir/wgs_reads/*
	do
		if  [ "$(cat $f | wc -l)" -ne 4000000 ]
		then
			echo "$f"
			echo "not enough reads"
			rm $f
		else
			file_name="${f%%.*}"
			bn="${file_name##*/}"
			profile="$dir/profiles/$bn.profile"
			if [ -f $profile ]
			then
				continue
			else
				motus profile -s $f -o $profile -t 25 -C precision -l 45 -g 1
			fi
		fi
	done
	#conda deactivate
done
