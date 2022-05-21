# WGSUniFrac-reproducibles
### Pre-requisites

* conda or miniconda
* python 3.6 (or higher)

### Installation & dependencies

```
git clone https://github.com/KoslickiLab/WGSUniFrac-reproducibles.git
cd WGSUniFrac-reproducibles
conda env create -f data/wgsunifrac.yml
conda env create -f data/qiime2.yml
```

### 1.1 On taxonomic data converted from phylogenetic data

1. Acquire data 

```
mkdir results
mkdir results/exp1
bash get_data_exp1.sh
```

2. Generate raw data (16S biom tables and WGS profiles).

```
conda activate wgsunifrac
python generate_rawdata_exp1.py -o results/exp1 -dm data/sorted_distance_complete.txt -mf data/otu_with_valid_taxid.txt
```

3. 


