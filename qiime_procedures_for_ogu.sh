#!/bin/bash
input_dir=$1
qiime tools import --type FeatureTable[Frequency] --input-path "$input_dir"/ogu.biom --output-path "$input_dir"/ogu.qza
qiime diversity core-metrics-phylogenetic --i-phylogeny tree.qza --i-table "$input_dir"/ogu.qza --p-sampling-depth 1000 --m-metadata-file "$input_dir"/metadata.tsv --output-dir "$input_dir"/ogu_out
qiime tools export --input-path "$input_dir"/ogu_out/weighted_unifrac_pcoa_results.qza --output-path "$input_dir"/exported
qiime tools export --input-path "$input_dir"/ogu_out/weighted_unifrac_distance_matrix.qza --output-path "$input_dir"/exported
qiime tools export --input-path "$input_dir"/ogu_out/weighted_unifrac_emperor.qzv --output-path "$input_dir"/exported


