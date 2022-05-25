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
conda activate qiime2
bash qiime_procedures_GTDB.sh results/GTDB/exp1/testDissimilarity 
bash qiime_procedures_GTDB.sh results/GTDB/exp1/testRange 
conda activate wgsunifrac
python get_combined_dataframe.py -d results/GTDB/exp1/testDissimilarity -a -1 -s 'results/GTDB/exp1/testDissimilarity_combined_df.txt' -t GTDB1
python get_combined_dataframe.py -d results/GTDB/exp1/testRange -a -1 -s 'results/GTDB/exp1/testRange_combined_df.txt' -t GTDB1
python generate_plot.py -d data/results/GTDB/exp1/testDissimilarity -a -1 -s "results/GTDB/exp1/testDissimilarity_combined_df.txt" 
python generate_plot.py -d data/results/GTDB/exp1/testRange -a -1 -s "results/GTDB/exp1/testRange_combined_df.txt" 
```

#### 1.3.2 GTDB taxonomy vs. NCBI taxonomy

### 4. OGU vs. WGSUniFrac

This section assumes Grinder is properly installed and working in the conda environment **wgsunifrac**.

1. Data preparation.

```
source get_data_OGU.sh 
```

2. Generate reads. Computationally heavy loaded. Modify do_get_wgs_libraries_ogu.sh according to the capacity of your machine.

```
for dir in results/OGU/testRange/*; do ./do_get_wgs_libraries_ogu.sh $dir; done 
for dir in results/OGU/testDissimilarity/*; do ./do_get_wgs_libraries_ogu.sh $dir; done
```

3. Get alignment files.

```
for dir in results/OGU/testRange/*; do ./do_get_alignment_file.sh $dir; done
for dir in results/OGU/testDissimilarity/*; do ./do_get_alignment_file.sh $dir; done
```

4. Woltka + Qiime beta diversity analysis (OGU procedures)

```
for dir in results/OGU/testRange/*; do woltka classify -i "$dir"/bt2out -m data/nucl2g.txt -o "$dir"/ogu.biom; done
for dir in results/OGU/testDissimilarity/*; do woltka classify -i "$dir"/bt2out -m data/nucl2g.txt -o "$dir"/ogu.biom; done
for dir in results/OGU/testRange/*; do ./mk_metadata.sh $dir; done #get metadata file
for dir in results/OGU/testDissimilarity/*; do ./mk_metadata.sh $dir; done #get metadata file
for dir in results/OGU/testRange/*; do ./qiime_procedures_for_ogu.sh $dir; done 
for dir in results/OGU/testDissimilarity/*; do ./qiime_procedures_for_ogu.sh $dir; done 
```

5. Get WGSUniFrac profiles.

```
for dir in results/OGU/testRange/*; do ./get_wgs_profiles.sh $dir; done #get profile
for dir in results/OGU/testDissimilarity/*; do ./get_wgs_profiles.sh $dir; done #get profile
```

6. Get plots.

```
python get_dataframe_ogu_vs_wgsunifrac.py -d results/OGU/testRange -s results/ogu_wgs_dataframe_env2_range.txt 
python get_dataframe_ogu_vs_wgsunifrac.py -d results/OGU/testDissimilarity -s results/ogu_wgs_dataframe_env2_dissimilarity.txt
python get_ogu_wgsunifrac_plot.py -f results/OGU/ogu_wgs_dataframe_env2_range.txt -x range -s "results/ogu_vs_wgsunifrac_range_lineplot.png"
python get_ogu_wgsunifrac_plot.py -f results/OGU/ogu_wgs_dataframe_env2_dissimilarity.txt -x range -s "results/ogu_vs_wgsunifrac_dissimilarity_lineplot.png"
```


