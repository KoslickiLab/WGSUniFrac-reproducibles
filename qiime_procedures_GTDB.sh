#!/bin/bash
input_dir=$1
qiime tools import --input-path data/bac120_r202.tree --output-path data/bac120_r202_tree.qza --type 'Phylogeny[Rooted]'
tree='data/bac120_r202_tree.qza'
#convert to .biom file
for exp_dir in $input_dir/*
do
	biom_file="$exp_dir/distance_matrix.biom"
	biom convert -i "$exp_dir/distance_matrix.txt" -o $biom_file --table-type="OTU table" --to-hdf5
	qiime tools import --input-path $biom_file --type 'FeatureTable[Frequency]' --input-format BIOMV210Format --output-path "$exp_dir/feature-table.qza"
	qiime diversity beta-phylogenetic --i-table "$exp_dir/feature-table.qza" --i-phylogeny $tree --p-metric weighted_unifrac --o-distance-matrix "$exp_dir/dmatrix.qza"
	qiime tools export --input-path "$exp_dir/dmatrix.qza" --output-path $exp_dir
done
