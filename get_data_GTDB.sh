#!/usr/bin/bash

#GTDB tree
wget -P ./data https://data.gtdb.ecogenomic.org/releases/release202/202.0/bac120_r202.tree
wget -P ./data https://data.gtdb.ecogenomic.org/releases/release202/202.0/bac120_metadata_r202.tar.gz
tar -C ./data -xvzf ./data/bac120_metadata_r202.tar.gz
wget -P ./data https://data.gtdb.ecogenomic.org/releases/release202/202.0/bac120_taxonomy_r202.tsv
wget -P ./data https://zenodo.org/record/6570678/files/bac120_distance_matrix_exp2.txt
wget -P ./data https://zenodo.org/record/6570678/files/bac120_valid_distance_matrix.txt

