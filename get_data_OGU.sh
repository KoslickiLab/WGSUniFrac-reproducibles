#!/bin/bash
#get WOL data
mkdir OGU
wget -P ./data/OGU -r -np -nH --cut-dirs=2 ftp://ftp.microbio.me/pub/wol-20April2021/
wget -P ./data/OGU https://biocore.github.io/wol/data/genomes/metadata.tsv.xz
conda activate qiime2 #for bowtie2
#build database
mkdir -p data/OGU/bowtie2
xz -d -k data/OGU/genomes/concat.fna.xz
qiime tools import --input-path data/OGU/trees/tree.nwk --output-path OGU/trees/WOL_tree.qza --type 'Phylogeny[Rooted]'
shuf data/OGU/taxonomy/taxid.map | head -3000 > data/OGU/taxonomy/taxid_map_mini.txt #shuffle and get 3000 random 
cat data/OGU/taxonomy/taxid_map_mini.txtt | cut -f1 > data/OGU/wol_nodes_3000.txt #get node list
head -1 metadata.tsv > metadata_mini.tsv
while read -r line; do grep $line data/OGU/metadata.tsv >> data/OGU/metadata_mini.tsv; done < data/OGU/wol_nodes_3000.txt 
python make_down_list.py data/OGU/metadata_mini.tsv > download.list
bash batch_down.sh download.list
mkdir data/wol_genome
mv G* data/wol_genome
#some files may be corrupted. remove.
for file in data/wol_genomes/*gz
do 
	echo $file
	file_name="${file%%.*}"
        bn="${file_name##*/}"
	gunzip $file 
	if ! [[ -s $file ]] #file empty
	then
		rm $file
		sed -i '/$bn/d' data/OGU/taxonomy/taxid_map_mini.txt
	fi
done
#concacanate 
touch data/OGU/genomes/wol_genome_mini.fasta
for file in data/wol_genomes/*; do echo $file; cat $file >> data/OGU/genomes/wol_genome_mini.fasta; done
rm -r data/wol_genomes
#build
bowtie2-build --thread 8 data/OGU/genomes/wol_genome_mini.fasta data/OGU/bowtie2/wol_mini
#get distance matrix
conda activate wgsunifrac
python get_distance_matrix.py -f data/OGU/taxonomy/taxid_map_mini.txt -t data/OGU/trees/tree.nwk -c 0 -gtdb 0 -o data/OGU/wol_distance_matrix.txt &

#create data
mkdir -p results/OGU
python get_abundance_files_for_ogu.py -od results/OGU -df data/OGU/wol_distance_matrix.txt -s 25 -o 200 -e 2


