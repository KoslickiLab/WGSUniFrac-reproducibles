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
python get_combined_dataframe.py -d results/exp1/testRange -a -1 -s 'results/exp1/combined_dataframe_range.txt' -t exp1
python get_combined_dataframe.py -d results/exp1/testDissimilarity -a -1 -s 'results/exp1/combined_dataframe_dissimilarity.txt' -t exp1
```

5. Get boxplots from the dataframe.

```
python generate_plot.py -f results/exp1/combined_dataframe_range.txt -x range -y Silhouette -s 'results/exp1/range_vs_silhouette_boxplot.png' 
python generate_plot.py -f results/exp1/combined_dataframe_dissimilarity.txt -x dissimilarity -y Silhouette -s 'results/exp1/dissimilarity_vs_silhouette_boxplot.png'
```

### 1.2 Model-based branch length assignment

The exponent of the branch length function defined can be adjusted by adjusting the -a parameter. Below is an example of setting the branch length function exponent to 0, i.e. all the branches are assigned the same length. All plots were generated in this manner.

```
python get_combined_dataframe.py -d results/exp1/testRange -a 0 -s 'results/exp1/combined_dataframe_range_bf0.txt' -t exp1
python get_combined_dataframe.py -d results/exp1/testDissimilarity -a 0 -s 'results/exp1/combined_dataframe_dissimilarity_bf0.txt' -t exp1
python generate_plot.py -f results/exp1/combined_dataframe_range_bf0.txt -x range -y silhouette -s 'results/exp1/range_vs_silhouette_boxplot.png' 
python generate_plot.py -f results/exp1/combined_dataframe_dissimilarity_bf0.txt -x dissimilarity -y silhouette -s 'results/exp1/dissimilarity_vs_silhouette_boxplot.png'
```

### 1.3 Data-driven branch lengths assignment using GTDB

1. Acquire data

```
conda activate wgsunifrac
bash get_data_GTDB.sh 
```

#### 1.3.1 GTDB vs. reciprocal model branch lengths assignment

```
mkdir results/GTDB
mkdir results/GTDB/exp1
python get_GTDB_input_1.py -od results/GTDB/exp1
bash qiime_procedures_GTDB.sh results/GTDB/exp1/testDissimilarity 
bash qiime_procedures_GTDB.sh results/GTDB/exp1/testRange &
python get_combined_dataframe.py -d results/GTDB/exp1/testDissimilarity -a -1 -s 'results/GTDB/exp1/testDissimilarity_combined_df.txt' -t GTDB1
python get_combined_dataframe.py -d results/exp1/GTDB/testRange -a -1 -s 'results/GTDB/exp1/testRange_combined_df.txt' -t GTDB1
python generate_plot.py -d data/results/GTDB/exp1/testDissimilarity -a -1 -s "results/GTDB/exp1/testDissimilarity_combined_df.txt" 
python generate_plot.py -d data/results/GTDB/exp1/testRange -a -1 -s "results/GTDB/exp1/testRange_combined_df.txt" 
```

#### 1.3.2 GTDB taxonomy vs. NCBI taxonomy


