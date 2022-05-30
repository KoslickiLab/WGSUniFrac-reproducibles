#!/bin/bash
input_dir=$1
#structure of directories:
#input_dir/exp_dir/16s_reads
for exp_dir in $input_dir/* 
do {
	#bash make_manifest.sh $exp_dir/16s_reads
	nohup qiime tools import --type 'SampleData[SequencesWithQuality]' --input-path $exp_dir/16s_reads/manifest.txt --output-path $exp_dir/demux.qza --input-format SingleEndFastqManifestPhred33
	nohup qiime dada2 denoise-single --i-demultiplexed-seqs $exp_dir/demux.qza --p-trunc-len 0 --output-dir $exp_dir/output --verbose
	nohup qiime fragment-insertion sepp --i-representative-sequences $exp_dir/output/representative_sequences.qza --i-reference-database /data/wjw5274/data/sepp-refs-gg-13-8.qza --o-tree $exp_dir/gg_13_8_tree.qza --o-placements $exp_dir/tree_placement.qza --p-threads 50
	nohup qiime diversity beta-phylogenetic --i-table $exp_dir/output/table.qza --i-phylogeny $exp_dir/gg_13_8_tree.qza --p-metric weighted_unifrac --o-distance-matrix $exp_dir/dmatrix.qza
	nohup qiime tools export --input-path $exp_dir/dmatrix.qza --output-path $exp_dir/exported_distance_matrix	
} &
done
