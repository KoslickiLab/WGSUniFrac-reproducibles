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
### Additional pre-requisites (install as needed)

* Grinder: for simulating amplicon and WGS reads. 
  * https://github.com/zyxue/biogrinder
* mOTUs: for profiling WGS reads into taxonomic profiles.
  * https://github.com/motu-tool/mOTUs

* wget
  * `brew install wget` with Homebrew or install via other methods

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

3. Get 16s pairwise UniFrac matrices using Qiime.

```
conda activate qiime2
bash get_16s_distance_matrix_with_qiime_exp1.sh results/exp1/testRange
bash get_16s_distance_matrix_with_qiime_exp1.sh results/exp1/testDissimilarity
```

4. Get combined dataframe.

```
conda activate wgsunifrac
python get_combined_dataframe_exp1.py -d results/exp1/testRange -a -1 -s 'results/exp1/combined_dataframe_range.txt'
python get_combined_dataframe_exp1.py -d results/exp1/testDissimilarity -a -1 -s 'results/exp1/combined_dataframe_dissimilarity.txt'
```

5. Get boxplots from the dataframe.

```
python generate_plot_exp1.py -f reproducibles/exp1/combined_dataframe_range.txt -x range -y silhouette -s 'reproducibles/exp1/range_vs_silhouette_boxplot.png' 
python generate_plot_exp1.py -f reproducibles/exp1/combined_dataframe_dissimilarity.txt -x dissimilarity -y silhouette -s 'reproducibles/exp1/dissimilarity_vs_silhouette_boxplot.png'
```


