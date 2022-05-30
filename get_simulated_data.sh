#!/bin/bash
[ -f data/grinder_distance_matrix_match_primer.txt ] || wget -P ./data https://zenodo.org/record/6587816/files/grinder_distance_matrix_match_primer.txt
mkdir -p results/simulated
python get_grinder_input.py -s 25 -o 200 -od results/simulated -e 2 
python get_grinder_input.py -s 25 -o 200 -od results/simulated -e 5
[ -f data/99_otus.fasta ] || wget -P ./data https://zenodo.org/record/6588302/files/99_otus.fasta
[ -f data/wgs_genome2.fasta ] || wget -P ./data https://zenodo.org/record/6588302/files/wgs_genome2.fasta
#get 16s libraries
bash get_16s_libraries_batch.sh results/simulated/testRange-env2
bash get_16s_libraries_batch.sh results/simulated/testRange-env5
bash get_16s_libraries_batch.sh results/simulated/testDissimilarity-env2
bash get_16s_libraries_batch.sh results/simulated/testDissimilarity-env5
bash get_wgs_libraries_batch.sh results/simulated/testRange-env2
bash get_wgs_libraries_batch.sh results/simulated/testRange-env5
bash get_wgs_libraries_batch.sh results/simulated/testDissimilarity-env2
bash get_wgs_libraries_batch.sh results/simulated/testDissimilarity-env5
for dir in results/simulated/testRange-env2; do bash make_manifest.sh $dir; done
for dir in results/simulated/testRange-env5; do bash make_manifest.sh $dir; done
for dir in results/simulated/testDissimilarity-env2; do bash make_manifest.sh $dir; done
for dir in results/simulated/testDissimilarity-env5; do bash make_manifest.sh $dir; done

